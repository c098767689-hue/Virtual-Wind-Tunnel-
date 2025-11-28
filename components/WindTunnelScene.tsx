
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

// --- SDF MATH FUNCTIONS (High Precision) ---

// Box SDF
const sdBox = (p: THREE.Vector3, b: THREE.Vector3) => {
  const dX = Math.abs(p.x) - b.x;
  const dY = Math.abs(p.y) - b.y;
  const dZ = Math.abs(p.z) - b.z;
  const inside = Math.min(Math.max(dX, Math.max(dY, dZ)), 0.0);
  const outside = Math.sqrt(Math.max(dX, 0) ** 2 + Math.max(dY, 0) ** 2 + Math.max(dZ, 0) ** 2);
  return inside + outside;
};

// Smooth Minimum (blends shapes)
const smin = (a: number, b: number, k: number) => {
  const h = Math.max(k - Math.abs(a - b), 0.0) / k;
  return Math.min(a, b) - h * h * k * (1.0 / 4.0);
};

// Capped Cylinder SDF
const sdCappedCylinder = (p: THREE.Vector3, h: number, r: number) => {
  const dXZ = Math.sqrt(p.x * p.x + p.z * p.z) - r;
  const dY = Math.abs(p.y) - h;
  return Math.min(Math.max(dXZ, dY), 0.0) + Math.sqrt(Math.max(dXZ, 0) ** 2 + Math.max(dY, 0) ** 2);
};

const getSDF = (pInput: { x: number; y: number; z: number }, type: ObjectType) => {
  // Reusable Vector3 to avoid allocation in tight loops? 
  // In JS, allocating small objects in a loop is okay, but for SDF thousands of times per frame, 
  // we try to keep math raw or simple objects.
  const p = { x: pInput.x, y: pInput.y, z: pInput.z };

  if (type === 'car') {
    // Car Body (Bottom)
    // Box: 1.8w, 0.5h, 3.8d -> Half: 0.9, 0.25, 1.9
    // Position: y = -0.25
    const pBody = { x: p.x, y: p.y + 0.25, z: p.z };
    const d1 = sdBox(pBody as THREE.Vector3, { x: 0.9, y: 0.25, z: 1.9 } as THREE.Vector3);
    
    // Cabin (Top)
    // Box: 1.4w, 0.5h, 1.8d -> Half: 0.7, 0.25, 0.9
    // Position: y = 0.25, z = -0.2
    const pCabin = { x: p.x, y: p.y - 0.25, z: p.z + 0.2 };
    const d2 = sdBox(pCabin as THREE.Vector3, { x: 0.7, y: 0.25, z: 0.9 } as THREE.Vector3);
    
    // Smooth blending to simulate aerodynamic curve
    return smin(d1, d2, 0.15);
  } 
  else if (type === 'cybertruck') {
    // We need to match the ExtrudeGeometry logic.
    // The visual mesh is rotated Y by -90 deg relative to standard.
    // In our physics loop, we handle 'cybertruck' by swapping x/z to align with the mesh's length.
    // Here, we define the shape assuming P is already in local "Truck Space" (Length along Z).
    
    // 1. Main lower body box
    // Width 2.0 (x +/- 1), Height 0.6, Length 4.0 (z +/- 2)
    const boxLower = sdBox({ x: p.x, y: p.y - 0.3, z: p.z } as THREE.Vector3, { x: 0.95, y: 0.3, z: 2.0 } as THREE.Vector3);
    
    // 2. The Triangle Roof/Bed
    // We can approximate the angular shape by intersecting planes or a rotated box.
    // Or simpler: A central triangular prism logic.
    
    // Sloped Plane 1 (Windshield): Normal pointing forward-up
    // Normal ~ [0, 1, 1] normalized
    const nz1 = 0.6; const ny1 = 0.8; // Approximate slope
    const distPlaneFront = (p.z * nz1 + p.y * ny1) - 0.8; // Offset
    
    // Sloped Plane 2 (Bed): Normal pointing back-up
    const nz2 = -0.5; const ny2 = 0.8;
    const distPlaneBack = (p.z * nz2 + p.y * ny2) - 0.8;

    // Intersection of Box and Planes (Rough approximation of Cybertruck peak)
    const boxUpper = sdBox({ x: p.x, y: p.y - 0.6, z: p.z } as THREE.Vector3, { x: 0.95, y: 0.6, z: 2.0 } as THREE.Vector3);
    
    // The truck is the intersection of the upper box bounds and the two cutting planes
    // But actually, smax (intersection) of the box and the space *below* the planes.
    // Actually, simpler: Union of lower box and a "Tent".
    
    // Let's use a simplified approach that fits inside the visual mesh:
    // A box that gets cut.
    const rawBox = sdBox({ x: p.x, y: p.y - 0.4, z: p.z } as THREE.Vector3, { x: 0.9, y: 0.6, z: 1.9 } as THREE.Vector3);
    
    // Cut front (Windshield)
    // Z > -2. Plane eqn: y - z - 1 = 0
    const cutFront = (p.y - (p.z * 0.45) - 0.9);
    
    // Cut back (Bed cover)
    // Z < 2. Plane eqn: y + z - 1 = 0
    const cutBack = (p.y + (p.z * 0.35) - 0.9);
    
    // Intersection: Max(Box, CutFront, CutBack)
    return Math.max(rawBox, cutFront, cutBack);
  } 
  else if (type === 'wing') {
    // Wing is rotated in physics loop so Length is X, Profile is YZ.
    // We want a teardrop shape in YZ.
    
    // 1. Limit length (X axis)
    const dX = Math.abs(p.x) - 2.2;
    
    // 2. Airfoil Profile in YZ
    // An airfoil is rounder at front (positive Z in local space?) and sharp at back.
    // Let's assume Z is flow direction (chord).
    // Modify Y based on Z.
    
    // Standard NACA-ish look:
    // Scale Y by a factor that depends on Z
    // Let's just do a stretched cylinder that tapers.
    const zNorm = (p.z + 0.5) / 1.5; // Normalize roughly 0 to 1 along chord
    const thickness = 0.25 * (1.0 - zNorm * 0.8); // Thicker at front
    const dProfile = Math.sqrt(p.y * p.y * 4.0) - thickness; // Elliptical approx
    
    // Simple Flattened Cylinder fallback if the above is unstable
    const ry = p.y / 0.15; // Thinner
    const rz = p.z / 1.0;  // Longer
    const dCyl = Math.sqrt(ry * ry + rz * rz) - 1.0;
    
    // Intersection of infinite profile and length cap
    return Math.max(dCyl * 0.15, dX);
  } 
  else if (type === 'sphere') {
    return Math.sqrt(p.x * p.x + p.y * p.y + p.z * p.z) - 1.0;
  } 
  else if (type === 'cylinder') {
    // Vertical cylinder: r=0.8, h=2 (half-height 1.0)
    return sdCappedCylinder(p as THREE.Vector3, 1.0, 0.8);
  } 
  else if (type === 'cone') {
    // Vertical Cone
    // r1 = 0.8 (bottom), r2 = 0.0 (top), h = 2.0
    // Simplified exact cone SDF is complex, using capped cone approx
    const q = Math.sqrt(p.x * p.x + p.z * p.z);
    // radius at y: from 0.8 at y=-1 to 0 at y=1
    // lerp: r = 0.4 * (1 - y)
    const r = 0.4 * (1.0 - p.y);
    const dSlope = (q - r) * 0.9; // cos(angle) factor approx
    const dHeight = Math.abs(p.y) - 1.0;
    return Math.max(dSlope, dHeight);
  } 
  else if (type === 'torus') {
    // Torus on XZ plane
    const t = { x: 1.0, y: 0.3 }; // major, minor
    const q = { x: Math.sqrt(p.x * p.x + p.z * p.z) - t.x, y: p.y };
    return Math.sqrt(q.x * q.x + q.y * q.y) - t.y;
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
        const y = 0;

        let lx = x * cosR - z * sinR;
        let ly = y;
        let lz = x * sinR + z * cosR;

        // MATCHING TRANSFORM LOGIC
        if (objectType === 'cybertruck') {
          // Mesh is rotated -90 Y. 
          // Physics is local. 
          // If we pass X,Z to physics, we need to swap them to match the "Long Z" definition in SDF
          const tmp = lx; lx = -lz; lz = tmp;
        } else if (objectType === 'wing') {
          // Wing rotates around X for AoA
          const radX = -angle * Math.PI / 180;
          const cosX = Math.cos(radX);
          const sinX = Math.sin(radX);
          
          // Revert the World Y Rotation for Wing (it stays at 0 yaw)
          // Actually, in TestObject, Wing has `rotation.y = 0`.
          // So we should NOT apply the `rad` rotation to `lx/lz` above.
          // Reset to world coords:
          lx = x; 
          const wy = y; 
          const wz = z;
          
          // Apply Pitch (X-axis rotation)
          ly = wy * cosX - wz * sinX;
          lz = wy * sinX + wz * cosX;
        }

        const dist = getSDF({ x: lx, y: ly, z: lz }, objectType);
        
        // ... Coloring Logic (Same as before) ...
        let cr = 0, cg = 0, cb = 0;
        let alpha = 0;

        if (dist < 0) {
           cr = 0.1; cg = 0.2; cb = 0.3; alpha = 0.9; // Solid inside
        } else {
           const d = dist; 
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

  useFrame(() => {
    if (!linesRef.current) return;
    
    const geom = linesRef.current.geometry;
    const posAttr = geom.attributes.position as THREE.BufferAttribute;
    const colAttr = geom.attributes.color as THREE.BufferAttribute;
    const renderPositions = posAttr.array as Float32Array;
    const renderColors = colAttr.array as Float32Array;
    
    const physPositions = physicsPosRef.current;
    const velocities = velocitiesRef.current;
    
    const rad = -angle * Math.PI / 180;
    const cosR = Math.cos(rad);
    const sinR = Math.sin(rad);

    let totalPressure = 0;
    let dragSum = 0;
    const targetSpeed = windSpeed * 0.4; // Tuned for better visual flow speed
    const trailFactor = 0.12 + (windSpeed * 0.1); 

    // Pre-calc Wing rotation
    const radX = -angle * Math.PI / 180;
    const cosX = Math.cos(radX);
    const sinX = Math.sin(radX);

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
        // Wing is static in Y, rotates in X (Pitch)
        lx = px;
        ly = py * cosX - pz * sinX;
        lz = py * sinX + pz * cosX;
      } else {
        // Standard Y-Rotation (Yaw)
        lx = px * cosR - pz * sinR;
        ly = py;
        lz = px * sinR + pz * cosR;

        if (objectType === 'cybertruck') {
          // Fix Coordinate Swap for Cybertruck Extrusion
          const tmp = lx; lx = -lz; lz = tmp;
        }
      }

      const dist = getSDF({x: lx, y: ly, z: lz}, objectType);
      
      // IMPROVED: Larger Interaction Distance for smoother avoidance
      const interactionDist = 0.8; 
      let pressure = 0;

      if (dist < interactionDist) {
        // Calculate Normal (Gradient)
        // Sample slightly wider for smoother normals on angular shapes
        const e = 0.05; 
        const nx = getSDF({x: lx + e, y: ly, z: lz}, objectType) - getSDF({x: lx - e, y: ly, z: lz}, objectType);
        const ny = getSDF({x: lx, y: ly + e, z: lz}, objectType) - getSDF({x: lx, y: ly - e, z: lz}, objectType);
        const nz = getSDF({x: lx, y: ly, z: lz + e}, objectType) - getSDF({x: lx, y: ly, z: lz - e}, objectType);
        
        const len = Math.sqrt(nx*nx + ny*ny + nz*nz);
        // Safety check for zero gradient (far field or center)
        const safeLen = len > 0 ? len : 1;
        const nnx = nx / safeLen; 
        const nny = ny / safeLen; 
        const nnz = nz / safeLen;

        // Transform Normal back to World Space
        let wnx, wny, wnz;

        if (objectType === 'wing') {
           wnx = nnx;
           wny = nny * cosX + nnz * sinX;
           wnz = -nny * sinX + nnz * cosX;
        } else {
           let localNx = nnx; let localNz = nnz;
           if (objectType === 'cybertruck') {
             const tmp = localNx; localNx = localNz; localNz = -tmp;
           }
           wnx = localNx * cosR + localNz * sinR;
           wny = nny;
           wnz = -localNx * sinR + localNz * cosR;
        }

        const dot = vx * wnx + vy * wny + vz * wnz;

        if (dist < 0.02) {
          // --- COLLISION RESPONSE (Inside) ---
          
          // 1. Hard position correction (push out)
          const pushOut = (0.02 - dist) + 0.01; // Extra nudge
          px += wnx * pushOut;
          py += wny * pushOut;
          pz += wnz * pushOut;

          // 2. Velocity Reflection or Slide
          if (dot < 0) {
            // Remove velocity into the wall (Slide)
            vx -= dot * wnx;
            vy -= dot * wny;
            vz -= dot * wnz;
            
            // Add Friction
            vx *= 0.9; vy *= 0.9; vz *= 0.9;

            pressure = 1.0;
            dragSum += 2.0;
          }
        } else {
          // --- AVOIDANCE RESPONSE (Outside but close) ---
          
          // Strength increases as we get closer
          const proximity = 1.0 - (dist / interactionDist);
          const steerFactor = Math.pow(proximity, 2) * 2.0; // Non-linear increase

          if (dot < 0) {
            // Steer velocity to be tangent
            vx -= dot * wnx * steerFactor; 
            vy -= dot * wny * steerFactor; 
            vz -= dot * wnz * steerFactor;
            
            // Slight repulsion to prevent hitting
            const repulsion = 0.01 * proximity;
            vx += wnx * repulsion;
            vy += wny * repulsion;
            vz += wnz * repulsion;
            
            pressure = proximity * 0.5;
            dragSum += pressure * 0.1;
          }
        }
      }

      // Global Wind Acceleration
      // Smoothly accelerate back to target wind speed
      const windInertia = 0.03;
      vz += (targetSpeed - vz) * windInertia;
      
      // Damping lateral movement to stabilize flow
      vx *= 0.98;
      vy *= 0.98;

      px += vx; py += vy; pz += vz;

      // Bounds Reset
      let reset = false;
      if (pz > 8 || Math.abs(px) > 10 || Math.abs(py) > 10) {
        pz = -8 - Math.random() * 4; // Spawn further back variation
        px = (Math.random() - 0.5) * 6;
        py = (Math.random() - 0.5) * 4;
        vx = 0; vy = 0; vz = targetSpeed;
        pressure = 0;
        reset = true;
      }

      // Update Physics Arrays
      physPositions[pIdx] = px;
      physPositions[pIdx+1] = py;
      physPositions[pIdx+2] = pz;
      velocities[pIdx] = vx;
      velocities[pIdx+1] = vy;
      velocities[pIdx+2] = vz;

      // --- RENDER MAPPING ---
      
      const speed = Math.sqrt(vx*vx + vy*vy + vz*vz);
      const normalizedSpeed = Math.min(speed / (targetSpeed * 1.5), 1.0);

      // Render Colors
      let r, g, b;
      if (pressure > 0.1) {
        // Impact / Pressure color (White/Red)
        const t = Math.min(pressure, 1.0);
        r = 1.0; g = 1.0 - t * 0.5; b = 1.0 - t;
        totalPressure += pressure;
      } else {
        // Flow color (Blue/Cyan/Green based on speed)
        // Slow (Greenish) -> Fast (Blue/Cyan)
        r = 0.0;
        g = 0.2 + (1.0 - normalizedSpeed) * 0.5; 
        b = 0.5 + normalizedSpeed * 0.5;
      }

      // Update Geometry
      renderPositions[rIdx+3] = px;
      renderPositions[rIdx+4] = py;
      renderPositions[rIdx+5] = pz;

      if (reset) {
         renderPositions[rIdx] = px;
         renderPositions[rIdx+1] = py;
         renderPositions[rIdx+2] = pz;
      } else {
         // Dynamic trail length based on speed
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

    // Metrics Update (Throttled by Frame)
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
    shape.moveTo(0, 0);
    shape.lineTo(2, 0);
    shape.lineTo(2, 0.8);
    shape.lineTo(1, 1.3);
    shape.lineTo(-2, 0.7);
    shape.lineTo(-2, 0);
    const extrudeSettings = { depth: 1.8, bevelEnabled: false };
    const geo = new THREE.ExtrudeGeometry(shape, extrudeSettings);
    geo.center();
    geo.rotateY(-Math.PI / 2); // Visual Rotation adjustment
    return geo;
  }, []);

  const wingGeo = useMemo(() => {
    const geo = new THREE.CylinderGeometry(1, 1, 4.5, 64);
    geo.scale(1, 0.2, 1);
    geo.rotateZ(Math.PI / 2);
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
        if (meshRef.current) meshRef.current.rotation.y = rRad - Math.PI/2; 
    } else if (type === 'wing') {
        if (meshRef.current) {
            meshRef.current.rotation.y = 0; // Wing stays straight in Yaw
            meshRef.current.rotation.x = rRad; // Wing pitches in X
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
