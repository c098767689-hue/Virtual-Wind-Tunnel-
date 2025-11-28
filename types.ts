
export type ObjectType = 'car' | 'cybertruck' | 'sphere' | 'cylinder' | 'cone' | 'torus';

export type TurbulencePreset = 'smooth' | 'chaotic' | 'high';

export interface SimulationState {
  objectType: ObjectType;
  windSpeed: number; // 0.1 to 1.0
  angle: number; // -180 to 180 degrees
  turbulencePreset: TurbulencePreset;
}

export interface SimulationMetrics {
  drag: number;
  lift: number;
  pressure: number;
}
