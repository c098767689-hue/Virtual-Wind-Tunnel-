
import React, { useRef, useMemo, useEffect } from 'react';
import { Canvas, useFrame } from '@react-three/fiber';
import { OrbitControls, Environment } from '@react-three/drei';
import * as THREE from 'three';
import { ObjectType } from '../types';

// Fix for missing JSX types in current environment
declare global {
  namespace JSX {
    interface IntrinsicElements {
      group: any;
      mesh: any;
      boxGeometry: any;
      meshStandardMaterial: any;
      points: any;
      bufferGeometry: any;
      bufferAttribute: any;
      pointsMaterial: any;
      lineSegments: any;
      lineBasicMaterial: any;
      fog: any;
      ambientLight: any;
      directionalLight: any;
      spotLight: any;
      gridHelper: any;
      sphereGeometry: any;
      cylinderGeometry: any;
      extrudeGeometry: any;
      coneGeometry: any;
      torusGeometry: any;
      planeGeometry: any;
      pointLight: any;
    }
  }
}

// --- MATH HELPERS (Raw Performance) ---

// Box SDF
const sdBox = (px: number, py: number, pz: number, bx: number, by: number, bz: number) => {
  const dx = Math.abs(px) - bx;
  const dy = Math.abs(py) - by;
  const dz = Math.abs(pz) - bz;
  const inside = Math.min(Math.max(dx, Math.max(dy, dz)), 0.0);
  const outside = Math.sqrt(Math.max(dx, 0) ** 2 + Math.max(dy, 0) ** 2 + Math.max(dz, 0) ** 2);
  return inside + outside;
};

// Smooth Minimum (blends shapes)
const smin = (a: number, b: number, k: number) => {
  const h = Math.max(k - Math.abs(a - b), 0.0) / k;
  return Math.min(a, b) - h * h * k * (1.0 / 4.0);
};

// Capped Cylinder SDF
const sdCappedCylinder = (px: number, py: number, pz: number, h: number, r: number) => {
  const dXZ = Math.sqrt(px * px + pz * pz) - r;
  const dY = Math.abs(py) - h;
  return Math.min(Math.max(dXZ, dY), 0.0) + Math.sqrt(Math.max(dXZ, 0) ** 2 + Math.max(dY, 0) ** 2);
};

// 2D Polygon SDF (Raw Math for performance)
// Verts: flat array [x0, y0, x1, y1, ...]
const sdPolygonRaw = (px: number, py: number, verts: number[]) => {
  let d = Infinity;
  let s = 1.0;
  const num = verts.length / 2;
  
  // Initial distance check (p - v0)^2
  const v0x = verts[0], v0y = verts[1];
  d = (px - v0x)**2 + (py - v0y)**2;

  for (let i = 0, j = num - 1; i < num; j = i, i++) {
    const vix = verts[i*2], viy = verts[i*2+1];
    const vjx = verts[j*2], vjy = verts[j*2+1];
    
    const ex = vjx - vix;
    const ey = vjy - viy;
    
    const wx = px - vix;
    const wy = py - viy;
    
    // clamp(dot(w,e)/dot(e,e), 0, 1)
    const dotWE = wx*ex + wy*ey;
    const dotEE = ex*ex + ey*ey;
    const t = Math.min(Math.max(dotWE / dotEE, 0.0), 1.0);
    
    const bx = wx - ex*t;
    const by = wy - ey*t;
    
    d = Math.min(d, bx*bx + by*by);
    
    // Winding number logic
    const c1 = py >= viy;
    const c2 = py < vjy;
    const c3 = ex * wy > ey * wx;
    
    if ((c1 && c2 && c3) || (!c1 && !c2 && !c3)) s *= -1.0;
  }
  
  return s * Math.sqrt(d);
};

// Ellipse Approximation (Gradient-based) for better flow over wings
const sdEllipseApprox = (y: number, z: number, ry: number, rz: number) => {
    // Simple scaling approximation is robust enough for this visual sim
    // d = (length(p/r) - 1) * min(r)
    // Avoid divide by zero
    const py = y / (ry > 0 ? ry : 0.01);
    const pz = z / (rz > 0 ? rz : 0.01);
    const len = Math.sqrt(py*py + pz*pz);
    return (len - 1.0) * Math.min(ry, rz);
};

// --- CURL NOISE GENERATOR (Trigonometric Approximation) ---
// Generates divergence-free noise field for realistic fluid swirls
const computeCurl = (x: number, y: number, z: number, time: number) => {
  const eps = 1e-4; // Finite difference epsilon

  // Potential function (Psi) - A composite of sin/cos waves
  // We need a vector potential (Psi_x, Psi_y, Psi_z)
  const potential = (px: number, py: number, pz: number) => {
    const t = time * 0.5;
    const scale1 = 1.5;
    const scale2 = 3.2;
    
    // Layer 1
    let valX = Math.sin(py * scale1 + t) + Math.cos(pz * scale1 - t);
    let valY = Math.sin(pz * scale1 + t) + Math.cos(px * scale1 - t);
    let valZ = Math.sin(px * scale1 + t) + Math.cos(py * scale1 - t);
    
    // Layer 2 (Detail)
    valX += 0.5 * Math.sin(pz * scale2 + t * 2.3);
    valY += 0.5 * Math.sin(px * scale2 + t * 2.3);
    valZ += 0.5 * Math.sin(py * scale2 + t * 2.3);

    return { x: valX, y: valY, z: valZ };
  };

  // Curl = (dPsi_z/dy - dPsi_y/dz, dPsi_x/dz - dPsi_z/dx, dPsi_y/dx - dPsi_x/dy)
  // Approximate derivatives
  const p0 = potential(x, y, z);
  const px_ = potential(x + eps, y, z);
  const py_ = potential(x, y + eps, z);
  const pz_ = potential(x, y, z + eps);

  const x_rate_y = (py_.x - p0.x) / eps;
  const x_rate_z = (pz_.x - p0.x) / eps;

  const y_rate_x = (px_.y - p0.y) / eps;
  const y_rate_z = (pz_.y - p0.y) / eps;

  const z_rate_x = (px_.z - p0.z) / eps;
  const z_rate_y = (py_.z - p0.z) / eps;

  return {
    x: z_rate_y - y_rate_z,
    y: x_rate_z - z_rate_x,
    z: y_rate_x - x_rate_y
  };
};

// --- GEOMETRY DATA ---
// Cybertruck Profile Vertices (YZ plane, Centered)
// Visual Height 1.3 -> Center 0.65.
// Vertices adjusted to be relative to center (0, 0.65, 0)
const CYBERTRUCK_PROFILE = [
    2.0, -0.65,  // Rear Bumper Bottom
    2.0, 0.15,   // Rear Bed Top
    1.0, 0.65,   // Peak
    -2.0, 0.05,  // Front Hood Top
    -2.0, -0.65  // Front Bumper Bottom
];

const getSDF = (p: { x: number; y: number; z: number }, type: ObjectType) => {
  if (type === 'car') {
    // Car Body (Bottom)
    const d1 = sdBox(p.x, p.y + 0.25, p.z, 0.9, 0.25, 1.9);
    // Cabin (Top)
    const d2 = sdBox(p.x, p.y - 0.25, p.z + 0.2, 0.7, 0.25, 0.9);
    // Smooth blending
    return smin(d1, d2, 0.15);
  } 
  else if (type === 'cybertruck') {
    // Extruded Polygon along X (Width)
    // Profile in YZ plane (Length Z, Height Y)
    
    // 1. 2D Profile Distance
    // Note: p.z maps to Profile X (Length), p.y maps to Profile Y (Height)
    const dPoly = sdPolygonRaw(p.z, p.y, CYBERTRUCK_PROFILE);
    
    // 2. Extrusion width (X axis)
    // Width 1.8 -> Half 0.9
    const dExtrude = Math.abs(p.x) - 0.9;
    
    // Intersection (Max) of Profile and Width
    const m = Math.max(dPoly, dExtrude);
    return m;
  } 
  else if (type === 'wing') {
    // Wing along X axis. Length 4.5.
    // Cross section in YZ. Elliptical/Airfoil.
    // Thickness Y=0.2 (Radius), Chord Z=1.0 (Radius).
    
    // 1. Infinite Elliptical Cylinder
    const dProfile = sdEllipseApprox(p.y, p.z, 0.2, 1.0);
    
    // 2. Cap Length (X axis)
    const dLength = Math.abs(p.x) - 2.25;
    
    // Intersection
    return Math.max(dProfile, dLength);
  } 
  else if (type === 'sphere') {
    return Math.sqrt(p.x * p.x + p.y * p.y + p.z * p.z) - 1.0;
  } 
  else if (type === 'cylinder') {
    return sdCappedCylinder(p.x, p.y, p.z, 1.0, 0.8);
  } 
  else if (type === 'cone') {
    // Vertical Cone: r=0.8 bottom, 0 top. h=2.
    // Approx
    const q = Math.sqrt(p.x * p.x + p.z * p.z);
    const r = 0.4 * (1.0 - p.y);
    const dSlope = (q - r) * 0.9; 
    const dHeight = Math.abs(p.y) - 1.0;
    return Math.max(dSlope, dHeight);
  } 
  else if (type === 'torus') {
    const tx = 1.0; const ty = 0.3;
    const qx = Math.sqrt(p.x * p.x + p.z * p.z) - tx;
    return Math.sqrt(qx * qx + p.y * p.y) - ty;
  }
  return 10.0;
};

// --- COMPONENTS ---

const SDFSlice = ({ objectType, angle }: { objectType: ObjectType; angle: number }) => {
  const pointsRef = useRef<THREE.Points>(null);
  const width = 10;
  const depth = 10;
  const cols = 120;
  const rows = 120;
  
  const [positions, initialColors] = useMemo(() => {
    const pos = new Float32Array(cols * rows * 3);
    const col = new Float32Array(cols * rows * 3);
    let i = 0;
    for (let r = 0; r < rows; r++) {
      for (let c = 0; c < cols; c++) {
        const x = (c / (cols - 1)) * width - width / 2;
        const z = (r / (rows - 1)) * depth - depth / 2;
        pos[i * 3] = x; pos[i * 3 + 1] = 0; pos[i * 3 + 2] = z;
        col[i * 3] = 0; col[i * 3 + 1] = 0; col[i * 3 + 2] = 0;
        i++;
      }
    }
    return [pos, col];
  }, []);

  useFrame(() => {
    if (!pointsRef.current) return;
    const colAttr = pointsRef.current.geometry.attributes.color as THREE.BufferAttribute;
    const renderColors = colAttr.array as Float32Array;
    
    const rad = -angle * Math.PI / 180;
    const cosR = Math.cos(rad);
    const sinR = Math.sin(rad);

    let i = 0;
    for (let r = 0; r < rows; r++) {
      for (let c = 0; c < cols; c++) {
        const x = (c / (cols - 1)) * width - width / 2;
        const z = (r / (rows - 1)) * depth - depth / 2;
        
        let lx = x * cosR - z * sinR;
        let ly = 0;
        let lz = x * sinR + z * cosR;

        if (objectType === 'wing') {
          // Wing is static Y, pitches in X.
          // Revert World Yaw (since Wing stays 0 yaw visual)
          lx = x; 
          const wz = z;
          
          const radX = -angle * Math.PI / 180;
          const cosX = Math.cos(radX);
          const sinX = Math.sin(radX);
          
          ly = 0 * cosX - wz * sinX;
          lz = 0 * sinX + wz * cosX;
        } 
        
        const dist = getSDF({ x: lx, y: ly, z: lz }, objectType);
        
        // Coloring Logic
        let cr = 0, cg = 0, cb = 0;
        let alpha = 0;

        if (dist < 0) {
           cr = 0.1; cg = 0.2; cb = 0.3; alpha = 0.9;
        } else {
           const d = dist; 
           // Heatmap Gradient
           if (d < 0.5) { cr = 1.0; cg = d * 2.0; cb = 0.0; } 
           else if (d < 1.0) { cr = 1.0 - (d - 0.5) * 2.0; cg = 1.0; cb = 0.0; } 
           else if (d < 2.0) { cr = 0.0; cg = 1.0 - (d - 1.0) * 0.5; cb = (d - 1.0) * 0.5; } 
           else { cr = 0.0; cg = 0.5 - (d - 2.0) * 0.2; cb = 0.5 + (d - 2.0) * 0.2; }

           const interval = 0.25;
           const lineW = 0.03;
           const mod = d % interval;
           if (mod < lineW || mod > (interval - lineW)) {
             cr += 0.3; cg += 0.3; cb += 0.3; alpha = 0.8;
           } else {
             alpha = Math.max(0, 0.6 - d * 0.2);
           }
        }
        
        renderColors[i * 3] = cr * alpha;
        renderColors[i * 3 + 1] = cg * alpha;
        renderColors[i * 3 + 2] = cb * alpha;
        i++;
      }
    }
    colAttr.needsUpdate = true;
  });

  return (
    <points ref={pointsRef} position={[0, -0.05, 0]}>
      <bufferGeometry>
        <bufferAttribute attach="attributes-position" count={positions.length / 3} array={positions} itemSize={3} />
        <bufferAttribute attach="attributes-color" count={initialColors.length / 3} array={initialColors} itemSize={3} />
      </bufferGeometry>
      <pointsMaterial size={0.08} vertexColors transparent opacity={0.6} sizeAttenuation depthWrite={false} />
    </points>
  );
};


const PARTICLE_COUNT = 8000;

const Particles = ({ 
  objectType, 
  windSpeed, 
  angle, 
  onMetricsUpdate 
}: { 
  objectType: ObjectType; 
  windSpeed: number; 
  angle: number; 
  onMetricsUpdate: (m: { drag: number; lift: number; pressure: number }) => void 
}) => {
  const linesRef = useRef<THREE.LineSegments>(null);
  const velocitiesRef = useRef<Float32Array>(new Float32Array(PARTICLE_COUNT * 3));
  const physicsPosRef = useRef<Float32Array>(new Float32Array(PARTICLE_COUNT * 3));

  const [positions, colors] = useMemo(() => {
    const pos = new Float32Array(PARTICLE_COUNT * 2 * 3);
    const col = new Float32Array(PARTICLE_COUNT * 2 * 3);
    const vels = velocitiesRef.current;
    const physPos = physicsPosRef.current;
    
    for (let i = 0; i < PARTICLE_COUNT; i++) {
      const x = (Math.random() - 0.5) * 6;
      const y = (Math.random() - 0.5) * 4;
      const z = (Math.random() - 0.5) * 12;

      physPos[i * 3] = x; physPos[i * 3 + 1] = y; physPos[i * 3 + 2] = z;
      vels[i * 3 + 2] = 0.5;

      const rIdx = i * 6;
      pos[rIdx] = x; pos[rIdx+1] = y; pos[rIdx+2] = z;
      pos[rIdx+3] = x; pos[rIdx+4] = y; pos[rIdx+5] = z;
    }
    return [pos, col];
  }, []);

  useFrame((state) => {
    if (!linesRef.current) return;
    
    const geom = linesRef.current.geometry;
    const posAttr = geom.attributes.position as THREE.BufferAttribute;
    const colAttr = geom.attributes.color as THREE.BufferAttribute;
    const renderPositions = posAttr.array as Float32Array;
    const renderColors = colAttr.array as Float32Array;
    
    const physPositions = physicsPosRef.current;
    const velocities = velocitiesRef.current;
    
    const time = state.clock.elapsedTime;
    
    const rad = -angle * Math.PI / 180;
    const cosR = Math.cos(rad);
    const sinR = Math.sin(rad);
    
    // Wing specific rotation (Pitch X)
    const radX = -angle * Math.PI / 180;
    const cosX = Math.cos(radX);
    const sinX = Math.sin(radX);

    let totalPressure = 0;
    let dragSum = 0;
    const targetSpeed = windSpeed * 0.4;
    const trailFactor = 0.12 + (windSpeed * 0.1); 

    for (let i = 0; i < PARTICLE_COUNT; i++) {
      const pIdx = i * 3;
      const rIdx = i * 6;
      
      let px = physPositions[pIdx];
      let py = physPositions[pIdx+1];
      let pz = physPositions[pIdx+2];
      
      let vx = velocities[pIdx];
      let vy = velocities[pIdx+1];
      let vz = velocities[pIdx+2];

      // TRANSFORM TO LOCAL SDF SPACE
      let lx, ly, lz;

      if (objectType === 'wing') {
        lx = px;
        ly = py * cosX - pz * sinX;
        lz = py * sinX + pz * cosX;
      } else {
        lx = px * cosR - pz * sinR;
        ly = py;
        lz = px * sinR + pz * cosR;
      }

      const dist = getSDF({x: lx, y: ly, z: lz}, objectType);
      
      const interactionDist = 0.5; 
      let pressure = 0;
      let inBoundaryLayer = false;
      let isWake = false;

      if (dist < interactionDist) {
        // Gradient (Normal) Calculation
        const e = 0.02; // More precise epsilon
        const nx = getSDF({x: lx + e, y: ly, z: lz}, objectType) - getSDF({x: lx - e, y: ly, z: lz}, objectType);
        const ny = getSDF({x: lx, y: ly + e, z: lz}, objectType) - getSDF({x: lx, y: ly - e, z: lz}, objectType);
        const nz = getSDF({x: lx, y: ly, z: lz + e}, objectType) - getSDF({x: lx, y: ly, z: lz - e}, objectType);
        
        let len = nx*nx + ny*ny + nz*nz;
        if (len > 0) {
            len = Math.sqrt(len);
            const nnx = nx / len; 
            const nny = ny / len; 
            const nnz = nz / len;

            // Transform Normal back to World Space
            let wnx, wny, wnz;

            if (objectType === 'wing') {
               wnx = nnx;
               wny = nny * cosX + nnz * sinX;
               wnz = -nny * sinX + nnz * cosX;
            } else {
               wnx = nnx * cosR + nnz * sinR;
               wny = nny;
               wnz = -nnx * sinR + nnz * cosR;
            }

            const dot = vx * wnx + vy * wny + vz * wnz;

            if (dist < 0.05) {
              // --- HARD COLLISION / IMPACT ---
              let pushOut = 0;
              if (dist <= 0) {
                  // Penetration correction: Snap to surface + slight push
                  pushOut = -dist + 0.03;
              } else {
                  // Cushion
                  pushOut = (0.05 - dist) * 0.4;
              }
              
              px += wnx * pushOut;
              py += wny * pushOut;
              pz += wnz * pushOut;

              if (dot < 0) {
                // Bounce/Slide - Friction Logic
                // Kill normal velocity (No bounce, fluid sticks or slides)
                vx -= dot * wnx * 1.05; // 1.05 = slight bounce
                vy -= dot * wny * 1.05;
                vz -= dot * wnz * 1.05;
                
                // Friction (Simulate viscosity near surface)
                vx *= 0.85; 
                vy *= 0.85; 
                vz *= 0.85;

                pressure = 1.0;
                dragSum += 2.0;
              }
            } else {
              // --- AVOIDANCE / LAMINAR FLOW ---
              const proximity = 1.0 - (dist / interactionDist);
              inBoundaryLayer = true;

              if (dot < 0) {
                // Steer Tangent
                const steerStrength = Math.pow(proximity, 3) * 2.0;
                
                // Remove velocity towards object
                vx -= dot * wnx * steerStrength; 
                vy -= dot * wny * steerStrength; 
                vz -= dot * wnz * steerStrength;
                
                // Add soft repulsion to maintain volume
                const repulsion = 0.005 * proximity;
                px += wnx * repulsion;
                py += wny * repulsion;
                pz += wnz * repulsion;
                
                pressure = proximity * 0.6;
                dragSum += pressure * 0.2;
              }
            }
        }
      }
      
      // --- ADVANCED TURBULENCE (CURL NOISE) ---
      
      // Calculate Curl Noise Vector for this position
      // Scale coordinates for noise frequency
      const curl = computeCurl(px * 0.5, py * 0.5, pz * 0.5, time);
      
      // 1. Wake Detection (Behind object)
      // Approximate bounding box check based on object size
      let zWakeStart = 1.5;
      let yWakeWidth = 1.0;
      let xWakeWidth = 1.5;
      if (objectType === 'wing') { zWakeStart = 1.0; yWakeWidth = 0.3; }

      // Check if particle is "shadowed" by the object
      if (pz > zWakeStart && Math.abs(px) < xWakeWidth && Math.abs(py) < yWakeWidth) {
         isWake = true;
      }
      
      // 2. Apply Forces
      if (isWake) {
         // --- WAKE TURBULENCE ---
         // Strong swirling, reduced forward velocity (Drag)
         const wakeIntensity = 0.04 * windSpeed;
         vx += curl.x * wakeIntensity;
         vy += curl.y * wakeIntensity;
         vz += curl.z * wakeIntensity;
         
         // Drag suction
         vz -= 0.005 * windSpeed; 
      } else if (inBoundaryLayer) {
         // --- BOUNDARY LAYER ---
         // Small micro-turbulences near surface
         const blIntensity = 0.01 * windSpeed;
         vx += curl.x * blIntensity;
         vy += curl.y * blIntensity;
         
         // Viscous drag (slow down near surface)
         vx *= 0.96;
         vy *= 0.96;
         vz *= 0.96;
      } else {
         // --- FREE STREAM ---
         // Very subtle ambient air movement
         const ambientIntensity = 0.002 * windSpeed;
         vx += curl.x * ambientIntensity;
         vy += curl.y * ambientIntensity;
      }

      // 3. Global Wind & Inertia
      const windInertia = isWake ? 0.01 : 0.05; // Wake has less connection to free stream
      vz += (targetSpeed - vz) * windInertia;
      
      // Damping lateral movement to keep stream vaguely straight
      vx *= 0.99;
      vy *= 0.99;

      // Integration
      px += vx; py += vy; pz += vz;

      // Reset
      let reset = false;
      const boundsX = 8;
      const boundsY = 5;
      const boundsZ = 12;
      
      if (pz > boundsZ || Math.abs(px) > boundsX || Math.abs(py) > boundsY) {
        pz = -8 - Math.random() * 4;
        px = (Math.random() - 0.5) * 6;
        py = (Math.random() - 0.5) * 4;
        
        // Initial velocity with slight random variation
        vx = (Math.random() - 0.5) * 0.01; 
        vy = (Math.random() - 0.5) * 0.01; 
        vz = targetSpeed;
        
        pressure = 0;
        reset = true;
      }

      physPositions[pIdx] = px;
      physPositions[pIdx+1] = py;
      physPositions[pIdx+2] = pz;
      velocities[pIdx] = vx;
      velocities[pIdx+1] = vy;
      velocities[pIdx+2] = vz;

      // Render Update
      const speed = Math.sqrt(vx*vx + vy*vy + vz*vz);
      const normalizedSpeed = Math.min(speed / (targetSpeed * 1.5), 1.0);

      let r, g, b;
      if (pressure > 0.1) {
        const t = Math.min(pressure * 1.2, 1.0);
        r = 1.0; g = 1.0 - t * 0.5; b = 1.0 - t;
        totalPressure += pressure;
      } else {
        // Velocity Color Map
        // Fast (Laminar) = Cyan/Blue
        // Slow (Wake/Stagnation) = Deep Blue/Purple
        r = 0.0;
        g = 0.2 + (1.0 - normalizedSpeed) * 0.6; // More green if slow
        b = 0.5 + normalizedSpeed * 0.5;
        
        if (isWake) {
           // Tint wake particles slightly purple/darker
           r = 0.2; g *= 0.5; b *= 0.8;
        }
      }

      renderPositions[rIdx+3] = px;
      renderPositions[rIdx+4] = py;
      renderPositions[rIdx+5] = pz;

      if (reset) {
         renderPositions[rIdx] = px;
         renderPositions[rIdx+1] = py;
         renderPositions[rIdx+2] = pz;
      } else {
         const lag = trailFactor * 2.5;
         renderPositions[rIdx] = px - vx * lag;
         renderPositions[rIdx+1] = py - vy * lag;
         renderPositions[rIdx+2] = pz - vz * lag;
      }

      renderColors[rIdx+3] = r;
      renderColors[rIdx+4] = g;
      renderColors[rIdx+5] = b;

      const tailOpacity = Math.max(0, 0.5 * normalizedSpeed);
      renderColors[rIdx] = r * tailOpacity;
      renderColors[rIdx+1] = g * tailOpacity;
      renderColors[rIdx+2] = b * tailOpacity;
    }

    posAttr.needsUpdate = true;
    colAttr.needsUpdate = true;

    if (onMetricsUpdate) {
        const lift = Math.sin(angle * Math.PI / 180) * windSpeed * 2000;
        onMetricsUpdate({
            drag: (dragSum / (PARTICLE_COUNT * windSpeed)) * 10 + 0.1,
            lift: lift,
            pressure: totalPressure
        });
    }
  });

  return (
    <lineSegments ref={linesRef}>
      <bufferGeometry>
        <bufferAttribute attach="attributes-position" count={positions.length / 3} array={positions} itemSize={3} />
        <bufferAttribute attach="attributes-color" count={colors.length / 3} array={colors} itemSize={3} />
      </bufferGeometry>
      <lineBasicMaterial vertexColors transparent opacity={0.8} blending={THREE.AdditiveBlending} depthWrite={false} linewidth={1} />
    </lineSegments>
  );
};

const TestObject = ({ type, angle, showWireframe }: { type: ObjectType; angle: number; showWireframe: boolean }) => {
  const meshRef = useRef<THREE.Mesh>(null);
  const groupRef = useRef<THREE.Group>(null);

  const material = useMemo(() => new THREE.MeshStandardMaterial({
    color: 0x1e293b,
    roughness: 0.2,
    metalness: 0.9,
    side: THREE.DoubleSide
  }), []);

  useEffect(() => {
    material.wireframe = showWireframe;
  }, [showWireframe, material]);

  // Geometries
  const carGeo = useMemo(() => new THREE.BoxGeometry(1.8, 0.5, 3.8), []);
  const carCabinGeo = useMemo(() => new THREE.BoxGeometry(1.4, 0.5, 1.8), []);
  
  const cyberTruckGeo = useMemo(() => {
    const shape = new THREE.Shape();
    // Profile matches vertices used in SDF (shifted by +0.65 to be 0..1.3 range relative to bottom)
    // Vertices in SDF are Centered Y (-0.65 to 0.65). 
    // Shape uses coordinates relative to anchor. ExtrudeGeometry centers it later.
    shape.moveTo(0, 0);
    shape.lineTo(2, 0);
    shape.lineTo(2, 0.8);
    shape.lineTo(1, 1.3);
    shape.lineTo(-2, 0.7);
    shape.lineTo(-2, 0);
    const extrudeSettings = { depth: 1.8, bevelEnabled: false };
    const geo = new THREE.ExtrudeGeometry(shape, extrudeSettings);
    geo.center(); 
    geo.rotateY(-Math.PI / 2); // Align Width to X, Length to Z
    return geo;
  }, []);

  const wingGeo = useMemo(() => {
    // Correct Order: 
    // 1. Create Cylinder (Axis Y).
    // 2. Rotate Z 90 -> Axis X.
    // 3. Scale -> X (Length) stays 1x (4.5), Y (Thickness) becomes 0.2x, Z (Chord) stays 1x.
    const geo = new THREE.CylinderGeometry(1, 1, 4.5, 64);
    geo.rotateZ(Math.PI / 2);
    geo.scale(1, 0.2, 1);
    return geo;
  }, []);

  const sphereGeo = useMemo(() => new THREE.SphereGeometry(1, 48, 48), []);
  const cylinderGeo = useMemo(() => new THREE.CylinderGeometry(0.8, 0.8, 2, 48), []);
  const coneGeo = useMemo(() => new THREE.ConeGeometry(0.8, 2, 32), []); 
  const torusGeo = useMemo(() => new THREE.TorusGeometry(1, 0.3, 16, 48), []); 

  useFrame(() => {
    const rRad = -angle * Math.PI / 180;
    
    if (type === 'car') {
        if (groupRef.current) groupRef.current.rotation.y = rRad;
    } else if (type === 'cybertruck') {
        // Fix: Removed extra rotation. The geometry is already aligned to Width X, Length Z.
        // Just apply angle of attack.
        if (meshRef.current) meshRef.current.rotation.y = rRad; 
    } else if (type === 'wing') {
        if (meshRef.current) {
            meshRef.current.rotation.y = 0; 
            meshRef.current.rotation.x = rRad;
        }
    } else {
        if (meshRef.current) meshRef.current.rotation.y = rRad;
    }
  });

  if (type === 'car') {
    return (
      <group ref={groupRef}>
        <mesh geometry={carGeo} material={material} position={[0, -0.25, 0]} />
        <mesh geometry={carCabinGeo} material={material} position={[0, 0.25, -0.2]} />
      </group>
    );
  }

  let geo: THREE.BufferGeometry = sphereGeo;
  if (type === 'cybertruck') geo = cyberTruckGeo;
  if (type === 'wing') geo = wingGeo;
  if (type === 'cylinder') geo = cylinderGeo;
  if (type === 'cone') geo = coneGeo;
  if (type === 'torus') geo = torusGeo;

  return <mesh ref={meshRef} geometry={geo} material={material} />;
};

export const WindTunnelScene = ({
  objectType,
  windSpeed,
  angle,
  showWireframe,
  showSDF,
  onMetricsUpdate
}: {
  objectType: ObjectType;
  windSpeed: number;
  angle: number;
  showWireframe: boolean;
  showSDF: boolean;
  onMetricsUpdate: (m: { drag: number; lift: number; pressure: number }) => void;
}) => {
  return (
    <div className="absolute inset-0 z-0">
      <Canvas camera={{ position: [5, 3, 6], fov: 45 }}>
        <color attach="background" args={['#080808']} />
        <fog attach="fog" args={['#080808', 2, 25]} />
        
        {/* Advanced Lighting Environment */}
        <Environment preset="city" />
        
        {/* Overhead Strip Lights */}
        <group position={[0, 4, 0]}>
             <mesh position={[0, 0, -2]} rotation={[Math.PI/2, 0, 0]}>
                 <planeGeometry args={[8, 0.5]} />
                 <meshStandardMaterial color="#ffffff" emissive="#ffffff" emissiveIntensity={2} />
             </mesh>
             <mesh position={[0, 0, 2]} rotation={[Math.PI/2, 0, 0]}>
                 <planeGeometry args={[8, 0.5]} />
                 <meshStandardMaterial color="#ffffff" emissive="#ffffff" emissiveIntensity={2} />
             </mesh>
             <pointLight position={[0, -1, -2]} intensity={0.5} color="#ffffff" />
             <pointLight position={[0, -1, 2]} intensity={0.5} color="#ffffff" />
        </group>

        {/* Cinematic Main Light */}
        <spotLight 
          position={[-5, 5, 5]} 
          angle={0.4} 
          penumbra={0.5} 
          intensity={2} 
          color="#38bdf8" 
          castShadow 
        />
        <ambientLight intensity={0.1} />
        
        <gridHelper args={[20, 40, 0x333333, 0x111111]} position={[0, -1.5, 0]} />
        
        <TestObject type={objectType} angle={angle} showWireframe={showWireframe} />
        <Particles 
          objectType={objectType} 
          windSpeed={windSpeed} 
          angle={angle} 
          onMetricsUpdate={onMetricsUpdate} 
        />
        
        {showSDF && <SDFSlice objectType={objectType} angle={angle} />}
        
        <OrbitControls />
      </Canvas>
    </div>
  );
};
