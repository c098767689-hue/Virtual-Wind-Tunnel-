
import React, { useState, useRef, useCallback } from 'react';
import { WindTunnelScene } from './components/WindTunnelScene';
import { UIOverlay } from './components/UIOverlay';
import { ObjectType } from './types';

const App: React.FC = () => {
  const [objectType, setObjectType] = useState<ObjectType>('car');
  const [windSpeed, setWindSpeed] = useState<number>(0.5);
  const [angle, setAngle] = useState<number>(0);
  const [showWireframe, setShowWireframe] = useState<boolean>(false);
  const [showSDF, setShowSDF] = useState<boolean>(false);

  // We use a Ref to store metrics to avoid re-rendering the entire React tree
  // every frame (60fps) just to update some text numbers.
  const metricsRef = useRef({ drag: 0, lift: 0, pressure: 0 });

  const handleMetricsUpdate = useCallback((metrics: { drag: number; lift: number; pressure: number }) => {
    // Smooth the values slightly
    metricsRef.current.drag += (metrics.drag - metricsRef.current.drag) * 0.1;
    metricsRef.current.lift = metrics.lift; // Lift is calculated deterministically based on angle in this sim
    metricsRef.current.pressure += (metrics.pressure - metricsRef.current.pressure) * 0.1;
  }, []);

  return (
    <div className="w-full h-screen relative bg-black overflow-hidden">
      <WindTunnelScene 
        objectType={objectType}
        windSpeed={windSpeed}
        angle={angle}
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
