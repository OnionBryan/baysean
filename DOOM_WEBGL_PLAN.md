# DOOM-style WebGL JS Rebuild Plan (Three.js-first)

This is a practical blueprint to rebuild a classic **DOOM-like** shooter in JavaScript with **WebGL**, optimized for "buttery" feel (stable frametimes, low input latency, smooth camera and movement).

> Note: The original DOOM source is GPL-2.0. If you reuse code assets or logic directly from id Software repositories, keep licensing obligations in mind.

## 1) Stack decision: Three.js + raw math core

- **Renderer**: Three.js (WebGL2 backend) for rapid iteration and robust material/geometry pipeline.
- **Math core**: Your own matrix/vector transforms for simulation and camera control (deterministic and auditable).
- **Optional d3**: Use d3 only for development diagnostics (frametime histograms, heatmaps, profiling overlays), not runtime rendering.

Why this split:
- Three.js accelerates asset/view pipeline.
- You maintain explicit control over movement/physics/collision math.
- d3 is excellent for perf analytics, not game rendering.

## 2) "Butter" architecture (fixed-step simulation + interpolated render)

Use a fixed timestep simulation and decoupled render loop:

- Simulation: `dt = 1/120` (or 1/144) seconds fixed.
- Render: `requestAnimationFrame` variable rate.
- Interpolation: render between previous and current sim state by `alpha = accumulator / dt`.

This avoids variable-timestep instability and keeps controls consistent.

```js
const SIM_DT = 1 / 120;
let acc = 0;
let prev = performance.now() * 0.001;

function frame(nowMs) {
  const now = nowMs * 0.001;
  let frameDt = now - prev;
  prev = now;

  // clamp to avoid death spirals on tab return
  frameDt = Math.min(frameDt, 0.05);
  acc += frameDt;

  while (acc >= SIM_DT) {
    stepSimulation(SIM_DT);
    acc -= SIM_DT;
  }

  const alpha = acc / SIM_DT;
  renderInterpolated(alpha);
  requestAnimationFrame(frame);
}
```

## 3) Coordinate system and matrix rigor

Define conventions early and never violate:

- Right-handed world
- `+Y` up, `+X` right, `-Z` forward (Three.js style camera-forward is typically `-Z`)
- Column-major `mat4` (match Three.js internals)

### Core matrices

- World transform: `M_world = T * R * S`
- View transform: `V = inverse(M_camera)`
- Projection: `P = perspective(fov, aspect, near, far)`
- Clip transform: `p_clip = P * V * M_world * p_model`

For FPS camera:

- Store yaw/pitch scalars.
- Build camera rotation matrix from yaw/pitch.
- Derive forward/right vectors from rotation matrix columns (or rows, depending convention).

Example (conceptual):

```txt
R_yaw   = rotY(yaw)
R_pitch = rotX(pitch)
R_cam   = R_yaw * R_pitch

forward = normalize(R_cam * [0, 0, -1, 0])
right   = normalize(R_cam * [1, 0,  0, 0])
```

## 4) Movement model that feels responsive

Use acceleration-based movement with friction on the horizontal plane:

- Desired velocity from WASD in camera-local frame.
- Accelerate toward target velocity with capped acceleration.
- Apply ground friction when no input.
- Apply air-control (reduced lateral accel) when not grounded.

Discrete update:

```txt
v_{t+dt} = v_t + a * dt
x_{t+dt} = x_t + v_{t+dt} * dt
```

For stability:
- Use semi-implicit Euler (velocity update before position).
- Keep `dt` fixed.
- Clamp extreme velocities.

## 5) DOOM-style level representation in modern JS

Two options:

1. **Authentic route**: parse WAD and BSP nodes/segs/subsectors.
2. **Practical route**: custom JSON sectors + linedefs, then compile to runtime structures.

Recommended for speed: start practical, then add WAD import.

Data structures:
- `Sector`: floor/ceiling heights, light, texture refs.
- `Linedef`: endpoints, sidedefs, flags (solid/portal/trigger).
- Spatial partition: BSP or uniform grid for collision and visibility.

## 6) Collision and queries (matrix-aware)

Player can be approximated by a vertical capsule or cylinder.

Pipeline:
1. Integrate intended velocity.
2. Broadphase query nearby segments/sectors.
3. Narrowphase: segment-distance and penetration normal.
4. Resolve by projecting velocity onto collision tangent plane:

```txt
v' = v - n * dot(v, n)    // slide response
```

For stairs/steps:
- Test a "step-up" candidate by offsetting position up to `stepHeight`.
- If clear, accept stepped solution.

## 7) Rendering path for high FPS

- Start with unlit or lightly lit materials (avoid expensive PBR until needed).
- Texture atlases reduce state changes.
- Merge static geometry by sector where possible.
- Use frustum culling + portal/BSP visibility culling.
- Use instancing for repeated props.

Three.js performance settings:
- `powerPreference: "high-performance"`
- Disable shadows initially.
- Keep post-processing minimal; if used, use lightweight FXAA only.

## 8) Input latency minimization

- Pointer lock for mouse look.
- Read raw input events into a small state buffer each frame.
- Apply input at the start of sim ticks.
- Avoid smoothing filters that add lag; prefer deterministic acceleration curves.

## 9) Net-new project skeleton

```txt
/src
  /core
    math/vec3.js
    math/mat4.js
    timing/fixedLoop.js
  /game
    playerController.js
    collision.js
    levelRuntime.js
    wad/
  /render
    renderer.js
    camera.js
    materials.js
  /debug
    perfOverlayD3.js
main.js
```

## 10) Milestone plan

1. **Week 1**: fixed loop, camera, WASD + mouse look, empty room.
2. **Week 2**: sector/linedef format + collision slide.
3. **Week 3**: BSP/visibility + textured walls/floor/ceiling.
4. **Week 4**: enemies/projectiles, hit-scan, pickups.
5. **Week 5**: polish pass for frametime consistency + latency.

## 11) "Feel like butter" acceptance metrics

Track these continuously:

- 95th percentile frametime < 10 ms on target machine.
- 99th percentile frametime < 14 ms.
- Input-to-photon median < 50 ms.
- No simulation divergence at fixed seed over 10-minute replay.

Use d3 dev dashboards for:
- Frametime histogram
- Jank spikes over time
- Input event to sim-application delay

## 12) Minimal starter checklist

- [ ] Fixed-step simulation loop implemented.
- [ ] Camera matrix from yaw/pitch verified.
- [ ] Movement acceleration/friction tuned.
- [ ] Segment collision + slide working.
- [ ] One sample level loaded.
- [ ] 120 FPS stress test with telemetry.

---

If you want, next step I can generate a **ready-to-run Vite + Three.js starter** with:
- fixed-step loop,
- matrix-based FPS camera,
- collision-ready player controller,
- d3 perf overlay panel.
