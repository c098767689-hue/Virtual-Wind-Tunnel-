
import React, { useEffect, useRef } from 'react';
import { ObjectType } from '../types';

interface UIOverlayProps {
  objectType: ObjectType;
  setObjectType: (t: ObjectType) => void;
  windSpeed: number;
  setWindSpeed: (s: number) => void;
  angle: number;
  setAngle: (a: number) => void;
  showWireframe: boolean;
  setShowWireframe: (b: boolean) => void;
  showSDF: boolean;
  setShowSDF: (b: boolean) => void;
  metricsRef: React.MutableRefObject<{ drag: number; lift: number; pressure: number }>;
}

export const UIOverlay: React.FC<UIOverlayProps> = ({
  objectType,
  setObjectType,
  windSpeed,
  setWindSpeed,
  angle,
  setAngle,
  showWireframe,
  setShowWireframe,
  showSDF,
  setShowSDF,
  metricsRef
}) => {
  const cdRef = useRef<HTMLSpanElement>(null);
  const liftRef = useRef<HTMLDivElement>(null);
  const pressRef = useRef<HTMLDivElement>(null);
  const barRef = useRef<HTMLDivElement>(null);

  // High frequency update loop for metrics to avoid React re-renders on every frame
  useEffect(() => {
    let frameId: number;
    const update = () => {
      const { drag, lift, pressure } = metricsRef.current;
      
      // Update DOM directly
      if (cdRef.current) cdRef.current.innerText = drag.toFixed(3);
      if (liftRef.current) liftRef.current.innerText = Math.floor(lift).toString();
      if (pressRef.current) pressRef.current.innerText = (pressure * 0.01).toFixed(2);
      
      if (barRef.current) {
        const width = Math.min((drag / 1.0) * 100, 100);
        barRef.current.style.width = `${width}%`;
      }
      
      frameId = requestAnimationFrame(update);
    };
    update();
    return () => cancelAnimationFrame(frameId);
  }, [metricsRef]);

  const objBtnClass = (active: boolean) => 
    `w-full p-2 rounded-lg text-xs font-bold transition-all border ${
      active 
      ? 'bg-sky-500 text-white border-sky-400 shadow-[0_0_15px_rgba(14,165,233,0.4)]' 
      : 'bg-slate-800 text-slate-400 border-slate-700 hover:bg-slate-700'
    }`;

  return (
    <div className="absolute inset-0 pointer-events-none font-sans select-none">
      
      {/* Top Header */}
      <div className="absolute top-0 w-full p-6 z-10 flex justify-between items-start">
        <div className="pointer-events-auto">
          <h1 className="text-3xl font-extrabold tracking-tight text-white flex items-center gap-3">
            <div className="w-8 h-8 bg-sky-500 rounded flex items-center justify-center text-black text-lg font-bold">W</div>
            <span>AERO<span className="text-sky-500">LAB</span></span>
          </h1>
          <div className="flex items-center gap-2 mt-2">
            <span className="text-[10px] bg-slate-800 text-slate-300 px-1.5 py-0.5 rounded border border-slate-700">SDF ENGINE</span>
            <span className="text-[10px] bg-slate-800 text-slate-300 px-1.5 py-0.5 rounded border border-slate-700">R3F CANVAS</span>
          </div>
        </div>
        
        <div className="px-4 py-2 rounded-full flex items-center gap-3 bg-slate-900/60 backdrop-blur-md border border-white/10 shadow-lg pointer-events-auto">
          <div className="flex flex-col items-end">
            <span className="text-[10px] text-slate-400 font-bold uppercase tracking-wider">Physics Rate</span>
            <span className="text-xs font-mono text-green-400">60 Hz</span>
          </div>
          <div className="h-6 w-[1px] bg-slate-700"></div>
          <div className="w-2 h-2 rounded-full bg-sky-500 animate-pulse shadow-[0_0_10px_#0ea5e9]"></div>
        </div>
      </div>

      {/* Left Controls */}
      <div className="absolute left-6 top-32 w-72 flex flex-col gap-4 pointer-events-none">
        
        {/* Model Selector */}
        <div className="p-5 rounded-2xl bg-slate-900/60 backdrop-blur-md border border-white/10 shadow-2xl pointer-events-auto transition-transform hover:scale-[1.02] duration-300">
          <div className="flex justify-between items-center mb-4 pb-2 border-b border-white/10">
            <h3 className="text-xs font-bold text-white uppercase tracking-widest">Test Model</h3>
          </div>
          <div className="grid grid-cols-3 gap-2 mb-4">
            <button onClick={() => setObjectType('car')} className={objBtnClass(objectType === 'car')}>Sports Car</button>
            <button onClick={() => setObjectType('cybertruck')} className={objBtnClass(objectType === 'cybertruck')}>Truck</button>
            <button onClick={() => setObjectType('wing')} className={objBtnClass(objectType === 'wing')}>Airfoil</button>
            <button onClick={() => setObjectType('sphere')} className={objBtnClass(objectType === 'sphere')}>Sphere</button>
            <button onClick={() => setObjectType('cylinder')} className={objBtnClass(objectType === 'cylinder')}>Cylinder</button>
            <button onClick={() => setObjectType('cone')} className={objBtnClass(objectType === 'cone')}>Cone</button>
            <button onClick={() => setObjectType('torus')} className={objBtnClass(objectType === 'torus')}>Torus</button>
          </div>
          
          <div className="grid grid-cols-2 gap-2">
            <button 
                onClick={() => setShowWireframe(!showWireframe)}
                className={`w-full py-1.5 px-3 rounded border text-[10px] font-bold uppercase tracking-wider transition-all flex items-center justify-center gap-2 ${
                showWireframe 
                    ? 'bg-sky-500/20 border-sky-500 text-sky-400' 
                    : 'bg-slate-800 border-slate-700 text-slate-400 hover:bg-slate-700'
                }`}
            >
                <span className={`w-2 h-2 rounded-full ${showWireframe ? 'bg-sky-400 shadow-[0_0_8px_#38bdf8]' : 'bg-slate-600'}`}></span>
                Wireframe
            </button>
            
            <button 
                onClick={() => setShowSDF(!showSDF)}
                className={`w-full py-1.5 px-3 rounded border text-[10px] font-bold uppercase tracking-wider transition-all flex items-center justify-center gap-2 ${
                showSDF 
                    ? 'bg-emerald-500/20 border-emerald-500 text-emerald-400' 
                    : 'bg-slate-800 border-slate-700 text-slate-400 hover:bg-slate-700'
                }`}
            >
                <span className={`w-2 h-2 rounded-full ${showSDF ? 'bg-emerald-400 shadow-[0_0_8px_#34d399]' : 'bg-slate-600'}`}></span>
                SDF Slice
            </button>
          </div>
        </div>

        {/* Sliders */}
        <div className="p-5 rounded-2xl bg-slate-900/60 backdrop-blur-md border border-white/10 shadow-2xl pointer-events-auto transition-transform hover:scale-[1.02] duration-300">
          <div className="flex justify-between items-center mb-4 pb-2 border-b border-white/10">
            <h3 className="text-xs font-bold text-white uppercase tracking-widest">Tunnel Controls</h3>
          </div>
          
          <div className="space-y-5">
            {/* Speed */}
            <div>
              <div className="flex justify-between text-xs mb-2">
                <span className="text-slate-400">Wind Velocity</span>
                <span className="text-sky-400 font-bold font-mono">{Math.round(windSpeed * 360)} km/h</span>
              </div>
              <input 
                type="range" 
                min="0.1" max="1.0" step="0.05" 
                value={windSpeed}
                onChange={(e) => setWindSpeed(parseFloat(e.target.value))}
              />
              <div className="flex justify-between text-[9px] text-slate-600 mt-1">
                <span>Low</span>
                <span>High</span>
              </div>
            </div>

            {/* Rotation */}
            <div>
              <div className="flex justify-between text-xs mb-2">
                <span className="text-slate-400">Angle of Attack</span>
                <span className="text-sky-400 font-bold font-mono">{angle}°</span>
              </div>
              <input 
                type="range" 
                min="-30" max="30" step="1" 
                value={angle}
                onChange={(e) => setAngle(parseFloat(e.target.value))}
              />
            </div>
          </div>
        </div>
      </div>

      {/* Right Data */}
      <div className="absolute right-6 top-32 w-64 flex flex-col gap-4 pointer-events-none">
        <div className="p-5 rounded-2xl bg-slate-900/60 backdrop-blur-md border border-white/10 shadow-2xl pointer-events-auto">
          <div className="flex justify-between items-center mb-4 pb-2 border-b border-white/10">
            <h3 className="text-xs font-bold text-white uppercase tracking-widest">Aerodynamics</h3>
            <div className="w-2 h-2 rounded-full bg-red-500 animate-pulse"></div>
          </div>

          <div className="space-y-4">
            <div className="bg-black/40 rounded-lg p-3 border border-white/5">
              <div className="text-[10px] text-slate-500 uppercase mb-1">Drag Coefficient</div>
              <div className="flex items-baseline gap-2">
                <span ref={cdRef} className="text-2xl font-bold text-white font-mono">0.000</span>
                <span className="text-xs text-slate-500">Cd</span>
              </div>
              <div className="w-full bg-slate-800 h-1.5 mt-2 rounded-full overflow-hidden">
                <div ref={barRef} className="h-full bg-gradient-to-r from-green-500 to-red-500 w-0 transition-all duration-300"></div>
              </div>
            </div>

            <div className="grid grid-cols-2 gap-2">
              <div className="bg-black/40 rounded-lg p-3 border border-white/5">
                <div className="text-[10px] text-slate-500 uppercase mb-1">Lift</div>
                <div ref={liftRef} className="text-lg font-bold text-white font-mono">0</div>
                <div className="text-[9px] text-slate-600">Newtons</div>
              </div>
              <div className="bg-black/40 rounded-lg p-3 border border-white/5">
                <div className="text-[10px] text-slate-500 uppercase mb-1">Pressure</div>
                <div ref={pressRef} className="text-lg font-bold text-orange-400 font-mono">0.00</div>
                <div className="text-[9px] text-slate-600">kPa (Max)</div>
              </div>
            </div>
          </div>
        </div>
      </div>

      {/* Legend */}
      <div className="absolute bottom-8 right-8 px-5 py-3 rounded-full bg-slate-900/60 backdrop-blur-md border border-white/10 flex items-center gap-4">
        <div className="flex items-center gap-2">
          <span className="w-3 h-3 rounded-full bg-blue-500 opacity-60"></span>
          <span className="text-[10px] text-slate-300 uppercase">Laminar</span>
        </div>
        <div className="flex items-center gap-2">
          <span className="w-3 h-3 rounded-full bg-white shadow-[0_0_10px_white]"></span>
          <span className="text-[10px] text-slate-300 uppercase">Impact</span>
        </div>
        <div className="flex items-center gap-2">
          <span className="w-3 h-3 rounded-full bg-red-500 opacity-80"></span>
          <span className="text-[10px] text-slate-300 uppercase">Turbulence</span>
        </div>
      </div>

    </div>
  );
};
