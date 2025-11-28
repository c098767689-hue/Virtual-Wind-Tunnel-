
import React, { useRef, useMemo, useEffect, useState } from 'react';
import { Canvas, useFrame } from '@react-three/fiber';
import { OrbitControls } from '@react-three/drei';
import * as THREE from 'three';
import { ObjectType, TurbulencePreset } from '../types';

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
      instancedMesh: any;
      cylinderGeometry: any;
      fog: any;
      ambientLight: any;
      directionalLight: any;
      spotLight: any;
      gridHelper: any;
      sphereGeometry: any;
      extrudeGeometry: any;
      coneGeometry: any;
      torusGeometry: any;
      planeGeometry: any;
      pointLight: any;
      hemisphereLight: any;
    }
  }
}

const FLOOR_LEVEL = -2.0;

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

// Capped Cone SDF (Exact)
const sdCappedCone = (px: number, py: number, pz: number, h: number, r1: number, r2: number) => {
  const qx = Math.sqrt(px * px + pz * pz);
  const qy = py;
  
  const k1x = r2; 
  const k1y = h;
  const k2x = r2 - r1; 
  const k2y = 2.0 * h;
  
  const cax = qx - Math.min(qx, (qy < 0.0) ? r1 : r2);
  const cay = Math.abs(qy) - h;
  
  const k2Len2 = k2x*k2x + k2y*k2y;
  const dotVal = (k1x - qx) * k2x + (k1y - qy) * k2y;
  const t = Math.min(Math.max(dotVal / k2Len2, 0.0), 1.0);
  
  const cbx = qx - k1x + k2x * t;
  const cby = qy - k1y + k2y * t;
  
  const s = (cbx < 0.0 && cay < 0.0) ? -1.0 : 1.0;
  
  return s * Math.sqrt(Math.min(cax*cax + cay*cay, cbx*cbx + cby*cby));
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

// --- CURL NOISE GENERATOR (Trigonometric Approximation) ---
// Generates divergence-free noise field for realistic fluid swirls
const computeCurl = (x: number, y: number, z: number, time: number, scaleFactor: number = 1.0) => {
  const eps = 1e-4; // Finite difference epsilon

  // Potential function (Psi) - A composite of sin/cos waves
  const potential = (px: number, py: number, pz: number) => {
    const t = time * 0.5;
    const scale1 = 1.5 * scaleFactor;
    const scale2 = 3.2 * scaleFactor;
    
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
const CYBERTRUCK_PROFILE = [
    2.0, -0.65,  // Rear Bumper Bottom
    2.0, 0.15,   // Rear Bed Top
    1.0, 0.65,   // Peak
    -2.0, 0.05,  // Front Hood Top
    -2.0, -0.65  // Front Bumper Bottom
];

const getSDF = (p: { x: number; y: number; z: number }, type: ObjectType) => {
  let dObject = 10.0;

  if (type === 'car') {
    // Car Body (Bottom)
    const d1 = sdBox(p.x, p.y + 0.25, p.z, 0.9, 0.25, 1.9);
    // Cabin (Top)
    const d2 = sdBox(p.x, p.y - 0.25, p.z + 0.2, 0.7, 0.25, 0.9);
    // Smooth blending
    dObject = smin(d1, d2, 0.15);
  } 
  else if (type === 'cybertruck') {
    // 1. 2D Profile Distance
    const dPoly = sdPolygonRaw(p.z, p.y, CYBERTRUCK_PROFILE);
    // 2. Extrusion width (X axis)
    const dExtrude = Math.abs(p.x) - 0.9;
    dObject = Math.max(dPoly, dExtrude);
  } 
  else if (type === 'sphere') {
    dObject = Math.sqrt(p.x * p.x + p.y * p.y + p.z * p.z) - 1.0;
  } 
  else if (type === 'cylinder') {
    dObject = sdCappedCylinder(p.x, p.y, p.z, 1.0, 0.8);
  } 
  else if (type === 'cone') {
    // Vertical Cone: Half-Height 1.0, Bottom Radius 0.8, Top Radius 0.0
    dObject = sdCappedCone(p.x, p.y, p.z, 1.0, 0.8, 0.0);
  } 
  else if (type === 'torus') {
    // Torus in XY plane (facing Z flow)
    const q = Math.sqrt(p.x * p.x + p.y * p.y) - 1.0;
    dObject = Math.sqrt(q * q + p.z * p.z) - 0.3;
  }

  // --- ADD FLOOR SDF ---
  // Simple plane at y = FLOOR_LEVEL
  const dFloor = p.y - FLOOR_LEVEL;
  
  // Return the union (min) of the object and the floor
  return Math.min(dObject, dFloor);
};

// --- COMPONENTS ---

const SDFSlice = ({ objectType, angle }: { objectType: ObjectType; angle: number }) => {
  const pointsRef = useRef<THREE.Points>(null);
  const width = 10;
  const depth = 10;
  const cols = 150; 
  const rows = 150;
  
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

  useFrame((state) => {
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
        
        const dist = getSDF({ x: lx, y: ly, z: lz }, objectType);
        
        // --- VISUALIZATION LOGIC ---
        let cr = 0, cg = 0, cb = 0;
        let alpha = 0;

        if (dist < 0) {
           // Inside Object
           cr = 0.8; cg = 0.1; cb = 0.2; 
           alpha = 0.8;
        } 
        else if (dist < 0.05) {
           // Boundary Line
           cr = 1.0; cg = 1.0; cb = 1.0; 
           alpha = 1.0;
        }
        else {
           const d = dist;
           if (d < 0.5) {
               const t = d / 0.5;
               cr = 1.0; cg = t; cb = 0.0;
           } else if (d < 1.0) {
               const t = (d - 0.5) / 0.5;
               cr = 1.0 - t; cg = 1.0; cb = 0.0;
           } else if (d < 2.0) {
               const t = (d - 1.0) / 1.0;
               cr = 0.0; cg = 1.0 - t * 0.5; cb = t;
           } else {
               const t = Math.min(1.0, (d - 2.0) / 2.0);
               cr = 0.0; cg = 0.5 * (1 - t); cb = 1.0 - t * 0.7;
           }

           // Contour Lines
           const interval = 0.5;
           const lineW = 0.04;
           const mod = d % interval;
           
           if (mod < lineW || mod > (interval - lineW)) {
             cr = cr * 0.5 + 0.5; 
             cg = cg * 0.5 + 0.5; 
             cb = cb * 0.5 + 0.5; 
             alpha = 0.9;
           } else {
             alpha = Math.max(0.1, 0.7 - d * 0.2);
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
      <pointsMaterial size={0.06} vertexColors transparent opacity={1.0} sizeAttenuation depthWrite={false} />
    </points>
  );
};


const PARTICLE_COUNT = 20000;
const MAX_SPARKS = 300;
const TRAIL_LENGTH = 5; // Segments in the tail

const Particles = ({ 
  objectType, 
  windSpeed, 
  angle, 
  turbulencePreset,
  onMetricsUpdate 
}: { 
  objectType: ObjectType; 
  windSpeed: number; 
  angle: number; 
  turbulencePreset: TurbulencePreset;
  onMetricsUpdate: (m: { drag: number; lift: number; pressure: number }) => void 
}) => {
  const meshRef = useRef<THREE.InstancedMesh>(null);
  
  // Physics State
  const velocitiesRef = useRef<Float32Array>(new Float32Array(PARTICLE_COUNT * 3));
  const physicsPosRef = useRef<Float32Array>(new Float32Array(PARTICLE_COUNT * 3));
  
  // Trail State
  const trailsRef = useRef<THREE.LineSegments>(null);
  const trailPosRef = useRef<Float32Array>(new Float32Array(PARTICLE_COUNT * TRAIL_LENGTH * 3));
  const trailColorRef = useRef<Float32Array>(new Float32Array(PARTICLE_COUNT * TRAIL_LENGTH * 3));
  
  // Sparks Data
  const sparkPointsRef = useRef<THREE.Points>(null);
  const sparkPosRef = useRef<Float32Array>(new Float32Array(MAX_SPARKS * 3));
  const sparkVelRef = useRef<Float32Array>(new Float32Array(MAX_SPARKS * 3));
  const sparkLifeRef = useRef<Float32Array>(new Float32Array(MAX_SPARKS)); // 0 = dead, 1 = max life
  const activeSparkIdx = useRef(0);

  const dummy = useMemo(() => new THREE.Object3D(), []);

  // Geometry for particle head
  const particleGeo = useMemo(() => {
    const geo = new THREE.SphereGeometry(0.025, 8, 8);
    return geo;
  }, []);

  // Geometry for trails (LineSegments)
  const trailGeometry = useMemo(() => {
    const geo = new THREE.BufferGeometry();
    const positions = new Float32Array(PARTICLE_COUNT * (TRAIL_LENGTH - 1) * 2 * 3); // 2 verts per segment
    const colors = new Float32Array(PARTICLE_COUNT * (TRAIL_LENGTH - 1) * 2 * 3);
    
    const indices = [];
    for (let i = 0; i < PARTICLE_COUNT; i++) {
        const offset = i * TRAIL_LENGTH;
        for (let j = 0; j < TRAIL_LENGTH - 1; j++) {
            indices.push(offset + j, offset + j + 1);
        }
    }
    
    geo.setIndex(indices);
    geo.setAttribute('position', new THREE.BufferAttribute(trailPosRef.current, 3));
    geo.setAttribute('color', new THREE.BufferAttribute(trailColorRef.current, 3));
    return geo;
  }, []);

  useEffect(() => {
    const physPos = physicsPosRef.current;
    const vels = velocitiesRef.current;
    const trailPos = trailPosRef.current;
    
    for (let i = 0; i < PARTICLE_COUNT; i++) {
        const x = (Math.random() - 0.5) * 6;
        // Ensure initial Y is above floor
        const y = Math.max(FLOOR_LEVEL + 0.1, (Math.random() - 0.5) * 4);
        const z = (Math.random() - 0.5) * 12;

        physPos[i * 3] = x; 
        physPos[i * 3 + 1] = y; 
        physPos[i * 3 + 2] = z;
        
        vels[i * 3 + 2] = 0.5 + Math.random() * 0.1;
        
        // Init trail history
        for (let t = 0; t < TRAIL_LENGTH; t++) {
            const tIdx = (i * TRAIL_LENGTH + t) * 3;
            trailPos[tIdx] = x;
            trailPos[tIdx+1] = y;
            trailPos[tIdx+2] = z;
        }
    }
  }, []);

  useFrame((state) => {
    if (!meshRef.current || !trailsRef.current) return;
    
    const physPositions = physicsPosRef.current;
    const velocities = velocitiesRef.current;
    const trailPos = trailPosRef.current;
    const trailCol = trailColorRef.current;
    
    const time = state.clock.elapsedTime;
    
    const rad = -angle * Math.PI / 180;
    const cosR = Math.cos(rad);
    const sinR = Math.sin(rad);

    let maxPressureFrame = 0;
    let dragSum = 0;
    const targetSpeed = windSpeed * 0.4;
    
    // PDF Section 2.1: Reynolds Number approximation
    const approxRe = Math.max(1.0, windSpeed * 2000); 
    const boundaryLayerThickness = 1.0 / Math.pow(approxRe, 0.4); 

    // PDF Section 3: Compressibility & Thermodynamics
    // If windSpeed > 0.8 (Approx Mach 0.3+ in simulator scale), apply heating effects
    const isCompressible = windSpeed > 0.8;
    const machHeatFactor = isCompressible ? (windSpeed - 0.8) * 5.0 : 0.0;

    let noiseScale = 1.0;
    let wakeMultiplier = 1.0;
    let ambientMultiplier = 1.0;
    let shedFreqMult = 1.0;

    switch (turbulencePreset) {
      case 'smooth':
        noiseScale = 0.4; // Reduced from 0.5 for smoother flow
        wakeMultiplier = 0.5; // Reduced from 0.8
        ambientMultiplier = 0.2; // Reduced from 0.5
        shedFreqMult = 0.5; // Reduced from 0.8
        break;
      case 'chaotic':
        noiseScale = 1.4; // Increased from 1.2
        wakeMultiplier = 3.5; // Increased from 2.5
        ambientMultiplier = 0.8; // Reduced slightly to focus on wake
        shedFreqMult = 2.0; // Increased from 1.5
        break;
      case 'high':
        noiseScale = 2.0; // Increased from 1.8
        wakeMultiplier = 2.0; // Increased from 1.5
        ambientMultiplier = 4.0; // Increased from 3.0
        shedFreqMult = 1.5; // Increased from 1.2
        break;
    }

    // --- SPARK SPAWN HELPER ---
    const spawnSpark = (x: number, y: number, z: number, nx: number, ny: number, nz: number) => {
        const idx = activeSparkIdx.current;
        activeSparkIdx.current = (activeSparkIdx.current + 1) % MAX_SPARKS;
        
        const pos = sparkPosRef.current;
        const vel = sparkVelRef.current;
        const life = sparkLifeRef.current;
        
        pos[idx * 3] = x;
        pos[idx * 3 + 1] = y;
        pos[idx * 3 + 2] = z;
        
        // Reflect velocity + random scatter
        const speed = 0.2 + Math.random() * 0.3;
        vel[idx * 3] = nx * speed + (Math.random() - 0.5) * 0.2;
        vel[idx * 3 + 1] = ny * speed + (Math.random() - 0.5) * 0.2;
        vel[idx * 3 + 2] = nz * speed + (Math.random() - 0.5) * 0.2;
        
        life[idx] = 1.0; // Full life
    };

    for (let i = 0; i < PARTICLE_COUNT; i++) {
      const pIdx = i * 3;
      
      let px = physPositions[pIdx];
      let py = physPositions[pIdx+1];
      let pz = physPositions[pIdx+2];
      
      let vx = velocities[pIdx];
      let vy = velocities[pIdx+1];
      let vz = velocities[pIdx+2];
      
      const speed = Math.sqrt(vx*vx + vy*vy + vz*vz);

      // Local Space Transform for SDF Check
      // Cars/Objects: Rotate around Y axis (Yaw)
      let lx = px * cosR - pz * sinR;
      let ly = py;
      let lz = px * sinR + pz * cosR;

      // Check object SDF + Floor SDF
      const dist = getSDF({x: lx, y: ly, z: lz}, objectType);
      
      const interactionDist = 0.6;
      let pressure = 0;
      let inBoundaryLayer = false;
      let isWake = false;

      // --- 1. PHYSICS INTERACTION ---
      if (dist < interactionDist) {
        // Gradient (Normal) Calculation
        const e = 0.02; 
        const nx = getSDF({x: lx + e, y: ly, z: lz}, objectType) - getSDF({x: lx - e, y: ly, z: lz}, objectType);
        const ny = getSDF({x: lx, y: ly + e, z: lz}, objectType) - getSDF({x: lx, y: ly - e, z: lz}, objectType);
        const nz = getSDF({x: lx, y: ly, z: lz + e}, objectType) - getSDF({x: lx, y: ly, z: lz - e}, objectType);
        
        let len = nx*nx + ny*ny + nz*nz;
        if (len > 0) {
            len = Math.sqrt(len);
            const nnx = nx / len; 
            const nny = ny / len; 
            const nnz = nz / len;

            // Transform Normal to World Space
            let wnx = nnx * cosR + nnz * sinR;
            let wny = nny;
            let wnz = -nnx * sinR + nnz * cosR;

            const dot = vx * wnx + vy * wny + vz * wnz;

            if (dist <= 0) {
                // Penetration - HARD CLAMP
                const penetration = Math.abs(dist);
                const pushOut = penetration + 0.01; // Increased safety margin
                
                // Immediate position correction to surface
                px += wnx * pushOut;
                py += wny * pushOut;
                pz += wnz * pushOut;
                
                // Update local transform var for next frame consistency
                lx += nnx * pushOut;
                ly += nny * pushOut;
                lz += nnz * pushOut;

                if (dot < 0) {
                    const impactIntensity = Math.abs(dot) / (speed + 0.001); 
                    const restitution = 0.6 + 0.35 * (1.0 - Math.min(1.0, impactIntensity)); 

                    const bounceFactor = 1.0 + restitution;
                    vx -= bounceFactor * dot * wnx;
                    vy -= bounceFactor * dot * wny;
                    vz -= bounceFactor * dot * wnz;

                    const scatterStrength = 0.15 * speed;
                    vx += (Math.random() - 0.5) * scatterStrength;
                    vy += (Math.random() - 0.5) * scatterStrength;
                    vz += (Math.random() - 0.5) * scatterStrength;

                    const impulse = Math.abs(dot);
                    pressure = impulse * 50.0;
                    dragSum += impulse;
                    
                    if (impulse > 0.05 && Math.random() < 0.2) {
                        spawnSpark(px, py, pz, wnx, wny, wnz);
                    }
                }
            } else if (dist < boundaryLayerThickness) {
                inBoundaryLayer = true;

                if (dot < 0) {
                    vx -= dot * wnx;
                    vy -= dot * wny;
                    vz -= dot * wnz;
                    
                    const impulse = Math.abs(dot);
                    pressure = impulse * 20.0 * (1.0 - dist/boundaryLayerThickness);
                    dragSum += impulse * 0.2; 
                }

                // PDF Section 1.2: Advanced Boundary Layer Theory (Log Law approximation)
                // y+ logic approximated by distance / thickness
                // Friction increases logarithmically as we get closer to 0 (No-Slip Condition)
                const yPlus = Math.max(0.001, dist / boundaryLayerThickness);
                const friction = Math.min(0.99, 0.1 + 0.15 * Math.log(yPlus * 100)); // Log-law approx
                
                vx *= friction;
                vy *= friction;
                vz *= friction;
            } else {
                // Pressure Field (Steering before impact)
                // If moving towards object (dot < 0), apply repulsive force
                if (dot < -0.1) {
                    const steerFactor = 0.02 * (1.0 - dist/interactionDist);
                    vx += wnx * steerFactor;
                    vy += wny * steerFactor;
                    vz += wnz * steerFactor;
                }
            }
        }
      }
      
      if (pressure > maxPressureFrame) {
          maxPressureFrame = pressure;
      }

      // --- 2. TURBULENCE & WAKE (RANS/LES Concepts) ---
      // PDF Section 1: Turbulence Modeling. 
      // We modify curl noise based on Shear Stress (velocity gradients).
      // Approximated here: Turbulence is higher where speed changes rapidly (boundary layer edge or wake).
      
      let localTurbulence = noiseScale;
      if (inBoundaryLayer) localTurbulence *= 2.0; // Higher shear in boundary layer

      const curl = computeCurl(px * 0.5, py * 0.5, pz * 0.5, time, localTurbulence);
      
      // Define Wake Zone
      let zWakeStart = 1.5;
      let baseW = 1.0; 
      let baseH = 0.8; 

      switch(objectType) {
        case 'car':
        case 'cybertruck':
            zWakeStart = 1.8; baseW = 1.0; baseH = 0.8;
            break;
        case 'sphere':
        case 'cylinder':
        case 'torus':
        case 'cone':
            zWakeStart = 0.8; baseW = 0.9; baseH = 0.9;
            break;
      }

      const distZ = pz - zWakeStart;
      
      if (distZ > 0) {
          const expansion = 1.0 + distZ * 0.1; // Reduced expansion
          const wakeW = baseW * expansion;
          const wakeH = baseH * expansion;
          
          if (Math.abs(px) < wakeW && Math.abs(py) < wakeH) {
             isWake = true;
          }
      }
      
      if (isWake) {
         const wakeIntensity = (0.1 + distZ * 0.05) * windSpeed * wakeMultiplier;
         
         vx += curl.x * wakeIntensity;
         vy += curl.y * wakeIntensity;
         vz += curl.z * wakeIntensity;
         
         // Von Karman Vortex Shedding (LES concept - Large Eddy Simulation)
         const shedFreq = 6.0 * windSpeed * shedFreqMult;
         const shedPhase = pz * 1.0 - time * shedFreq; 
         const vortexStrength = 0.2 * windSpeed * wakeMultiplier;
         
         vx += Math.sin(shedPhase) * vortexStrength;
         vy += Math.cos(shedPhase * 0.8) * vortexStrength * 0.6;
         
         if (distZ < 2.0) {
             const suctionFade = 1.0 - (distZ / 2.0);
             // Reduced wake suction significantly to prevent accumulation
             vz -= 0.01 * windSpeed * wakeMultiplier * suctionFade; 
         }
      } else if (!inBoundaryLayer) {
         const ambientIntensity = 0.002 * windSpeed * ambientMultiplier;
         vx += curl.x * ambientIntensity;
         vy += curl.y * ambientIntensity;
      }

      // --- 3. INTEGRATION ---
      const windInertia = isWake ? 0.02 : 0.05;
      vz += (targetSpeed - vz) * windInertia;
      vx *= 0.99; 
      vy *= 0.99;

      px += vx; py += vy; pz += vz;

      // --- 4. BOUNDS & RESET ---
      const boundsX = 8;
      const boundsY = 5;
      const boundsZ = 12;
      // Stagnation check
      const isStagnant = speed < 0.05 && pz > -2 && pz < 4 && Math.abs(px) < 2; 
      
      let didReset = false;
      if (pz > boundsZ || Math.abs(px) > boundsX || Math.abs(py) > boundsY || isStagnant) {
        pz = -8 - Math.random() * 4;
        px = (Math.random() - 0.5) * 6;
        // Reset above floor
        py = Math.max(FLOOR_LEVEL + 0.1, (Math.random() - 0.5) * 4);
        vx = (Math.random() - 0.5) * 0.01; 
        vy = (Math.random() - 0.5) * 0.01; 
        vz = targetSpeed;
        pressure = 0;
        didReset = true;
      }

      physPositions[pIdx] = px;
      physPositions[pIdx+1] = py;
      physPositions[pIdx+2] = pz;
      velocities[pIdx] = vx;
      velocities[pIdx+1] = vy;
      velocities[pIdx+2] = vz;

      // --- TRAIL UPDATE ---
      const trailOffset = i * TRAIL_LENGTH * 3;
      
      if (didReset) {
        for (let t = 0; t < TRAIL_LENGTH; t++) {
            const idx = trailOffset + t * 3;
            trailPos[idx] = px;
            trailPos[idx+1] = py;
            trailPos[idx+2] = pz;
        }
      } else {
        for (let t = TRAIL_LENGTH - 1; t > 0; t--) {
            const curr = trailOffset + t * 3;
            const prev = trailOffset + (t - 1) * 3;
            trailPos[curr] = trailPos[prev];
            trailPos[curr+1] = trailPos[prev+1];
            trailPos[curr+2] = trailPos[prev+2];
        }
        trailPos[trailOffset] = px;
        trailPos[trailOffset+1] = py;
        trailPos[trailOffset+2] = pz;
      }

      // --- RENDER UPDATE ---
      dummy.position.set(px, py, pz);
      dummy.lookAt(px + vx, py + vy, pz + vz);
      
      // No stretch for spheres
      dummy.scale.set(1, 1, 1);
      
      dummy.updateMatrix();
      meshRef.current.setMatrixAt(i, dummy.matrix);
      
      // Determine Color
      const normalizedSpeed = Math.min(speed / (targetSpeed * 1.5), 1.0);
      let r = 0, g = 0.5, b = 1.0;
      
      if (pressure > 0.1) {
          const t = Math.min(pressure * 0.2, 1);
          r = 1.0; g = 1.0 - t * 0.5; b = 1.0 - t;
      } else if (isCompressible) {
          // Thermodynamic heating visualization (Mach > 0.3)
          // Speed > 0.8 -> Color shifts to Orange/Red
          r = Math.min(1.0, machHeatFactor * 0.5);
          g = 0.5 + machHeatFactor * 0.2;
          b = 1.0 - machHeatFactor;
      } else if (isWake) {
          r = 0.6; g = 0.1; b = 0.8; 
      } else {
          r = 0.0; 
          g = normalizedSpeed * 0.8;
          b = 0.8 + normalizedSpeed * 0.2;
      }
      meshRef.current.setColorAt(i, new THREE.Color(r, g, b));

      for (let t = 0; t < TRAIL_LENGTH; t++) {
         const idx = trailOffset + t * 3;
         const opacity = 1.0 - (t / TRAIL_LENGTH);
         trailCol[idx] = r * opacity;
         trailCol[idx+1] = g * opacity;
         trailCol[idx+2] = b * opacity;
      }
    }

    meshRef.current.instanceMatrix.needsUpdate = true;
    if (meshRef.current.instanceColor) meshRef.current.instanceColor.needsUpdate = true;
    
    if (trailsRef.current) {
        trailsRef.current.geometry.attributes.position.needsUpdate = true;
        trailsRef.current.geometry.attributes.color.needsUpdate = true;
    }

    if (sparkPointsRef.current) {
        const sPos = sparkPosRef.current;
        const sVel = sparkVelRef.current;
        const sLife = sparkLifeRef.current;
        const sGeo = sparkPointsRef.current.geometry;
        const sColor = sGeo.attributes.color as THREE.BufferAttribute;
        const sP = sGeo.attributes.position as THREE.BufferAttribute;
        
        for (let i = 0; i < MAX_SPARKS; i++) {
            if (sLife[i] > 0) {
                sPos[i*3] += sVel[i*3];
                sPos[i*3+1] += sVel[i*3+1];
                sPos[i*3+2] += sVel[i*3+2];
                sVel[i*3+1] -= 0.01;
                sVel[i*3] *= 0.95;
                sVel[i*3+2] *= 0.95;
                sLife[i] -= 0.05;
            } else {
                sPos[i*3] = 0; sPos[i*3+1] = -100; sPos[i*3+2] = 0;
            }
            const l = Math.max(0, sLife[i]);
            sColor.setXYZ(i, 1.0, l * 0.8, l * 0.2);
        }
        sColor.needsUpdate = true;
        sP.needsUpdate = true;
    }

    if (onMetricsUpdate) {
        const angleRad = angle * Math.PI / 180;
        const Cl = 2 * Math.PI * Math.sin(angleRad); 
        const liftForce = Cl * (windSpeed * windSpeed) * 3000;
        const vSq = Math.max(0.01, windSpeed * windSpeed);
        const cd = (dragSum / PARTICLE_COUNT / vSq) * 1.2; 

        onMetricsUpdate({
            drag: cd,
            lift: liftForce,
            pressure: maxPressureFrame
        });
    }
  });

  return (
    <group>
        <instancedMesh ref={meshRef} args={[particleGeo, undefined, PARTICLE_COUNT]}>
            <meshStandardMaterial 
                color="#ffffff" 
                transparent 
                opacity={0.5} 
                roughness={0.4} 
                metalness={0.8} 
                emissive="#000044"
                emissiveIntensity={0.5}
                blending={THREE.AdditiveBlending}
                depthWrite={false}
            />
        </instancedMesh>
        
        <lineSegments ref={trailsRef} geometry={trailGeometry}>
             <lineBasicMaterial 
                vertexColors 
                transparent 
                opacity={0.3} 
                blending={THREE.AdditiveBlending} 
                depthWrite={false} 
            />
        </lineSegments>
        
        <points ref={sparkPointsRef}>
            <bufferGeometry>
                <bufferAttribute 
                    attach="attributes-position" 
                    count={MAX_SPARKS} 
                    array={sparkPosRef.current} 
                    itemSize={3} 
                />
                <bufferAttribute 
                    attach="attributes-color" 
                    count={MAX_SPARKS} 
                    array={new Float32Array(MAX_SPARKS * 3).fill(1)} 
                    itemSize={3} 
                />
            </bufferGeometry>
            <pointsMaterial 
                size={0.08} 
                vertexColors 
                transparent 
                opacity={1} 
                blending={THREE.AdditiveBlending} 
                depthWrite={false}
            />
        </points>
    </group>
  );
};

const TestObject = ({ type, angle, showWireframe }: { type: ObjectType; angle: number; showWireframe: boolean }) => {
  const meshRef = useRef<THREE.Mesh>(null);
  const groupRef = useRef<THREE.Group>(null);
  const [hovered, setHover] = useState(false);

  const material = useMemo(() => new THREE.MeshStandardMaterial({
    color: 0x1e293b,
    roughness: 0.2,
    metalness: 0.5,
    side: THREE.DoubleSide
  }), []);

  useEffect(() => {
    material.wireframe = showWireframe;
  }, [showWireframe, material]);

  useEffect(() => {
    document.body.style.cursor = hovered ? 'pointer' : 'auto';
    if (hovered) {
        material.emissive.setHex(0x38bdf8);
        material.emissiveIntensity = 0.2;
        material.color.setHex(0x2d3e54); 
    } else {
        material.emissive.setHex(0x000000);
        material.emissiveIntensity = 0.0;
        material.color.setHex(0x1e293b); 
    }
  }, [hovered, material]);

  const carGeo = useMemo(() => new THREE.BoxGeometry(1.8, 0.5, 3.8), []);
  const carCabinGeo = useMemo(() => new THREE.BoxGeometry(1.4, 0.5, 1.8), []);
  
  const cyberTruckGeo = useMemo(() => {
    const shape = new THREE.Shape();
    shape.moveTo(0, 0);
    shape.lineTo(2, 0);
    shape.lineTo(2, 0.8);
    shape.lineTo(1, 1.3);
    shape.lineTo(-2, 0.7);
    shape.lineTo(-2, 0);
    const extrudeSettings = { depth: 1.8, bevelEnabled: false };
    const geo = new THREE.ExtrudeGeometry(shape, extrudeSettings);
    geo.center(); 
    geo.rotateY(-Math.PI / 2);
    return geo;
  }, []);

  const sphereGeo = useMemo(() => new THREE.SphereGeometry(1, 48, 48), []);
  const cylinderGeo = useMemo(() => new THREE.CylinderGeometry(0.8, 0.8, 2, 48), []);
  const coneGeo = useMemo(() => new THREE.ConeGeometry(0.8, 2, 32), []); 
  const torusGeo = useMemo(() => {
    const geo = new THREE.TorusGeometry(1, 0.3, 16, 48);
    return geo;
  }, []); 

  useFrame(() => {
    const rRad = -angle * Math.PI / 180;
    
    if (type === 'car') {
        if (groupRef.current) groupRef.current.rotation.y = rRad;
    } else if (type === 'cybertruck') {
        if (meshRef.current) meshRef.current.rotation.y = rRad; 
    } else {
        if (meshRef.current) meshRef.current.rotation.y = rRad;
    }
  });

  if (type === 'car') {
    return (
      <group 
        ref={groupRef}
        onPointerOver={(e) => { e.stopPropagation(); setHover(true); }}
        onPointerOut={(e) => { e.stopPropagation(); setHover(false); }}
      >
        <mesh geometry={carGeo} material={material} position={[0, -0.25, 0]} />
        <mesh geometry={carCabinGeo} material={material} position={[0, 0.25, -0.2]} />
      </group>
    );
  }

  let geo: THREE.BufferGeometry = sphereGeo;
  if (type === 'cybertruck') geo = cyberTruckGeo;
  if (type === 'cylinder') geo = cylinderGeo;
  if (type === 'cone') geo = coneGeo;
  if (type === 'torus') geo = torusGeo;

  return (
    <mesh 
        ref={meshRef} 
        geometry={geo} 
        material={material} 
        onPointerOver={(e) => { e.stopPropagation(); setHover(true); }}
        onPointerOut={(e) => { e.stopPropagation(); setHover(false); }}
    />
  );
};

export const WindTunnelScene = ({
  objectType,
  windSpeed,
  angle,
  turbulencePreset,
  showWireframe,
  showSDF,
  onMetricsUpdate
}: {
  objectType: ObjectType;
  windSpeed: number;
  angle: number;
  turbulencePreset: TurbulencePreset;
  showWireframe: boolean;
  showSDF: boolean;
  onMetricsUpdate: (m: { drag: number; lift: number; pressure: number }) => void;
}) => {
  return (
    <div className="absolute inset-0 z-0">
      <Canvas camera={{ position: [5, 3, 6], fov: 45 }}>
        <color attach="background" args={['#080808']} />
        <fog attach="fog" args={['#080808', 2, 25]} />
        
        {/* Lights */}
        <hemisphereLight intensity={0.3} color="#ffffff" groundColor="#000000" />
        <spotLight 
          position={[-5, 5, 5]} 
          angle={0.4} 
          penumbra={0.5} 
          intensity={2} 
          color="#38bdf8" 
          castShadow 
        />
        <ambientLight intensity={0.1} />
        
        <gridHelper args={[20, 40, 0x333333, 0x111111]} position={[0, -2, 0]} />

        {/* Physical Floor Mesh */}
        <mesh rotation={[-Math.PI / 2, 0, 0]} position={[0, -2.01, 0]} receiveShadow>
            <planeGeometry args={[20, 40]} />
            <meshStandardMaterial color="#111111" roughness={0.8} metalness={0.2} />
        </mesh>
        
        <TestObject type={objectType} angle={angle} showWireframe={showWireframe} />
        <Particles 
          objectType={objectType} 
          windSpeed={windSpeed} 
          angle={angle} 
          turbulencePreset={turbulencePreset}
          onMetricsUpdate={onMetricsUpdate} 
        />
        
        {showSDF && <SDFSlice objectType={objectType} angle={angle} />}
        
        <OrbitControls />
      </Canvas>
    </div>
  );
};
