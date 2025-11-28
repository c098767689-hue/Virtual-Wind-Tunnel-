
import React, { useState, useRef, useCallback, useEffect } from 'react';
import { WindTunnelScene } from './components/WindTunnelScene';
import { UIOverlay } from './components/UIOverlay';
import { ObjectType, TurbulencePreset } from './types';

const STORAGE_KEY = 'aerolab_simulation_state_v1';

interface SavedState {
  objectType: ObjectType;
  windSpeed: number;
  angle: number;
  turbulencePreset: TurbulencePreset;
  showWireframe: boolean;
  showSDF: boolean;
}

const App: React.FC = () => {
  // Lazy initialization of state from localStorage
  // This runs only once on mount
  const [initialState] = useState<SavedState>(() => {
    const defaults: SavedState = {
      objectType: 'car',
      windSpeed: 0.5,
      angle: 0,
      turbulencePreset: 'smooth',
      showWireframe: false,
      showSDF: false
    };

    try {
      const saved = localStorage.getItem(STORAGE_KEY);
      if (saved) {
        const parsed = JSON.parse(saved);
        // Merge with defaults to ensure all fields exist
        return { ...defaults, ...parsed };
      }
    } catch (e) {
      console.warn('Failed to parse saved state:', e);
    }
    return defaults;
  });

  const [objectType, setObjectType] = useState<ObjectType>(initialState.objectType);
  const [windSpeed, setWindSpeed] = useState<number>(initialState.windSpeed);
  const [angle, setAngle] = useState<number>(initialState.angle);
  const [turbulencePreset, setTurbulencePreset] = useState<TurbulencePreset>(initialState.turbulencePreset);
  const [showWireframe, setShowWireframe] = useState<boolean>(initialState.showWireframe);
  const [showSDF, setShowSDF] = useState<boolean>(initialState.showSDF);

  // We use a Ref to store metrics to avoid re-rendering the entire React tree
  // every frame (60fps) just to update some text numbers.
  const metricsRef = useRef({ drag: 0, lift: 0, pressure: 0 });

  const handleMetricsUpdate = useCallback((metrics: { drag: number; lift: number; pressure: number }) => {
    // Smooth the values slightly
    metricsRef.current.drag += (metrics.drag - metricsRef.current.drag) * 0.1;
    metricsRef.current.lift = metrics.lift; // Lift is calculated deterministically based on angle in this sim
    metricsRef.current.pressure += (metrics.pressure - metricsRef.current.pressure) * 0.1;
  }, []);

  // Save state changes to localStorage
  useEffect(() => {
    const stateToSave: SavedState = {
      objectType,
      windSpeed,
      angle,
      turbulencePreset,
      showWireframe,
      showSDF
    };
    localStorage.setItem(STORAGE_KEY, JSON.stringify(stateToSave));
  }, [objectType, windSpeed, angle, turbulencePreset, showWireframe, showSDF]);

  return (
    <div className="w-full h-screen relative bg-black overflow-hidden">
      <WindTunnelScene 
        objectType={objectType}
        windSpeed={windSpeed}
        angle={angle}
        turbulencePreset={turbulencePreset}
        showWireframe={showWireframe}
        showSDF={showSDF}
        onMetricsUpdate={handleMetricsUpdate}
      />
      
      <UIOverlay 
        objectType={objectType}
        setObjectType={setObjectType}
        windSpeed={windSpeed}
        setWindSpeed={setWindSpeed}
        angle={angle}
        setAngle={setAngle}
        turbulencePreset={turbulencePreset}
        setTurbulencePreset={setTurbulencePreset}
        showWireframe={showWireframe}
        setShowWireframe={setShowWireframe}
        showSDF={showSDF}
        setShowSDF={setShowSDF}
        metricsRef={metricsRef}
      />
    </div>
  );
};

export default App;
