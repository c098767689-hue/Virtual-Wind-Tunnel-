export type ObjectType = 'car' | 'cybertruck' | 'wing' | 'sphere' | 'cylinder' | 'cone' | 'torus';

export interface SimulationState {
  objectType: ObjectType;
  windSpeed: number; // 0.1 to 1.0
  angle: number; // -30 to 30 degrees
}

export interface SimulationMetrics {
  drag: number;
  lift: number;
  pressure: number;
}