import React, { useEffect, useRef, useState } from 'react';
import * as $3DmolModule from '3dmol';

const $3Dmol = $3DmolModule.default || $3DmolModule;

const BUNDLE_URL = '/results/payload_unpacked/viz_bundle.json';
const UNPACKED_BASE = '/results/payload_unpacked/';
const N_RECEPTOR_VARIANTS = 6;
const RECEPTOR_FRAMES_PER_STEP = 8; // ligand frames per receptor conformer step (breathing rate)

// Ping-pong 0..5..0..5... so the receptor cycles smoothly through all 6 conformers and back,
// driven directly by the ligand's frame index (no separate animation timer to fall out of sync
// or misbehave on replay).
function receptorFrameForLigandFrame(frameIdx) {
    const period = 2 * (N_RECEPTOR_VARIANTS - 1);
    const pos = Math.floor(frameIdx / RECEPTOR_FRAMES_PER_STEP) % period;
    return pos <= N_RECEPTOR_VARIANTS - 1 ? pos : period - pos;
}

function xyzString(elements, positions) {
    const lines = [String(elements.length), 'frame'];
    for (let i = 0; i < elements.length; i++) {
        const [x, y, z] = positions[i];
        lines.push(`${elements[i]} ${x} ${y} ${z}`);
    }
    return lines.join('\n');
}

// Element reference for non-chemists: name + the exact Jmol/CPK color 3Dmol renders it in
// (3Dmol's default elemental coloring, used for every atom except ligand/pocket carbons which
// are recolored green/yellow below). Covers every element DiffSBDD samples during denoising.
const ELEMENT_LEGEND = [
    { sym: 'C', name: 'Carbon', color: '#909090' },
    { sym: 'N', name: 'Nitrogen', color: '#3050f8' },
    { sym: 'O', name: 'Oxygen', color: '#ff0d0d' },
    { sym: 'S', name: 'Sulfur', color: '#ffff30' },
    { sym: 'P', name: 'Phosphorus', color: '#ff8000' },
    { sym: 'F', name: 'Fluorine', color: '#90e050' },
    { sym: 'Cl', name: 'Chlorine', color: '#1ff01f' },
    { sym: 'Br', name: 'Bromine', color: '#a62929' },
    { sym: 'I', name: 'Iodine', color: '#940094' },
    { sym: 'B', name: 'Boron', color: '#ffb5b5' },
];
const ELEMENT_NAMES = Object.fromEntries(ELEMENT_LEGEND.map((e) => [e.sym, e.name]));

export default function TrajectoryViewer() {
    const containerRef = useRef(null);
    const viewerRef = useRef(null);
    const ligandModelRef = useRef(null);
    const receptorModelRef = useRef(null);
    const receptorFrameRef = useRef(-1);

    const [bundle, setBundle] = useState(null);
    const [error, setError] = useState(null);
    const [candidateIdx, setCandidateIdx] = useState(0);
    const [frameIdx, setFrameIdx] = useState(0);
    const [playing, setPlaying] = useState(true);
    const [speed, setSpeed] = useState(3);

    // Fetch the trajectory/candidate bundle once.
    useEffect(() => {
        fetch(BUNDLE_URL)
            .then((r) => {
                if (!r.ok) throw new Error(`Failed to load trajectory data (${r.status})`);
                return r.json();
            })
            .then((data) => {
                setBundle(data);
                const defaultIdx = data.candidates.findIndex((c) => c.default);
                setCandidateIdx(defaultIdx >= 0 ? defaultIdx : 0);
            })
            .catch((e) => setError(e.message));
    }, []);

    // Create the 3Dmol viewer and load the breathing receptor ensemble once the bundle is
    // known. The viewer is imperative (not React-reconciled), so it's built/torn down by hand
    // in this effect rather than via JSX.
    useEffect(() => {
        if (!bundle || !containerRef.current) return undefined;

        const viewer = $3Dmol.createViewer(containerRef.current, {
            backgroundColor: '#0b0f1a',
        });
        viewerRef.current = viewer;

        let cancelled = false;
        fetch(UNPACKED_BASE + bundle.meta.receptor_pdb)
            .then((r) => r.text())
            .then((pdbText) => {
                if (cancelled) return;
                const model = viewer.addModelsAsFrames(pdbText, 'pdb');
                receptorModelRef.current = model;
                model.setStyle({}, { cartoon: { color: '#5b7fb5', opacity: 0.85 } });

                // Highlight pocket-lining residues (precomputed by build_viz_bundle.py via
                // distance to the docking box center) so the cavity the ligand sits in reads
                // clearly, instead of just a uniform cartoon ribbon.
                const byChain = {};
                (bundle.meta.pocket_lining_residues || []).forEach(({ chain, resi }) => {
                    (byChain[chain] ||= []).push(resi);
                });
                Object.entries(byChain).forEach(([chain, resi]) => {
                    model.setStyle(
                        { chain, resi },
                        { cartoon: { color: '#f5a623' }, stick: { colorscheme: 'yellowCarbon', radius: 0.15 } }
                    );
                });

                viewer.zoomTo({ chain: Object.keys(byChain)[0], resi: Object.values(byChain)[0] });
                viewer.render();
            });

        return () => {
            cancelled = true;
            viewerRef.current = null;
            ligandModelRef.current = null;
            receptorModelRef.current = null;
            receptorFrameRef.current = -1;
            if (containerRef.current) containerRef.current.innerHTML = '';
        };
        // eslint-disable-next-line react-hooks/exhaustive-deps
    }, [bundle]);

    const candidate = bundle ? bundle.candidates[candidateIdx] : null;
    const nFrames = candidate ? candidate.frames.length : 0;

    // Swap the ligand model in place on every timeline tick, and step the receptor's displayed
    // conformer off the same frame index so the pocket visibly "breathes" in lockstep with the
    // denoising animation (including while scrubbing) rather than via a separate, independently
    // timed animation loop.
    //
    // Each denoising frame changes atom *types*, not just positions, so the ligand can't use
    // 3Dmol's addModelsAsFrames animation (it keeps atom identity fixed from frame 0) — instead
    // we remove/re-add a fresh model each tick, which makes 3Dmol re-perceive elements and
    // bonds, so bonds visibly form as the molecule resolves out of noise.
    useEffect(() => {
        const viewer = viewerRef.current;
        if (!viewer || !candidate) return;

        if (ligandModelRef.current) viewer.removeModel(ligandModelRef.current);
        const xyz = xyzString(candidate.elements[frameIdx], candidate.frames[frameIdx]);
        const model = viewer.addModel(xyz, 'xyz');
        model.setStyle({}, {
            stick: { radius: 0.18, colorscheme: 'greenCarbon' },
            sphere: { scale: 0.3, colorscheme: 'greenCarbon' },
        });
        ligandModelRef.current = model;

        // Hover-to-identify: atoms are re-added every tick (see comment above), so hoverable
        // state has to be re-armed on the current atom set each time too. Shows element name
        // for the ligand, element + residue for the receptor.
        viewer.setHoverable(
            {},
            true,
            (atom, v) => {
                if (atom.label) return;
                const name = ELEMENT_NAMES[atom.elem] || atom.elem;
                const text = atom.resn ? `${name} — ${atom.resn}${atom.resi}` : name;
                atom.label = v.addLabel(text, {
                    position: { x: atom.x, y: atom.y, z: atom.z },
                    backgroundColor: '#111827',
                    backgroundOpacity: 0.85,
                    fontColor: 'white',
                    fontSize: 12,
                    borderThickness: 0,
                });
            },
            (atom, v) => {
                if (atom.label) {
                    v.removeLabel(atom.label);
                    delete atom.label;
                }
            }
        );

        const receptorFrame = receptorFrameForLigandFrame(frameIdx);
        if (receptorModelRef.current && receptorFrame !== receptorFrameRef.current) {
            receptorFrameRef.current = receptorFrame;
            receptorModelRef.current.setFrame(receptorFrame).then(() => viewer.render());
        } else {
            viewer.render();
        }
    }, [candidate, frameIdx]);

    // Playback loop.
    useEffect(() => {
        if (!playing || !candidate) return undefined;
        const id = setInterval(() => {
            setFrameIdx((f) => (f + speed) % nFrames);
        }, 33);
        return () => clearInterval(id);
    }, [playing, speed, nFrames, candidate]);

    useEffect(() => {
        setFrameIdx(0);
    }, [candidateIdx]);

    if (error) {
        return (
            <div className="max-w-3xl w-full mx-auto mt-12 p-4 bg-red-50 border border-red-200 rounded-md">
                <p className="text-sm text-red-600">{error}</p>
                <p className="text-xs text-red-500 mt-2">
                    Run <code>venv/bin/python3 build_viz_bundle.py</code> to generate the trajectory bundle.
                </p>
            </div>
        );
    }

    if (!bundle || !candidate) {
        return <div className="text-center mt-12 text-gray-500">Loading denoising trajectories…</div>;
    }

    return (
        <div className="max-w-5xl w-full mx-auto mt-8 grid grid-cols-1 lg:grid-cols-3 gap-6">
            <div
                ref={containerRef}
                className="lg:col-span-2 bg-white border border-gray-200 rounded-lg shadow-sm overflow-hidden relative"
                style={{ height: 560 }}
            />

            <div className="space-y-4">
                <div className="bg-white p-4 border border-gray-200 shadow-sm rounded-lg">
                    <h3 className="text-sm font-medium text-gray-900 mb-2">Candidate</h3>
                    <select
                        className="w-full border border-gray-300 rounded-md px-2 py-1.5 text-sm"
                        value={candidateIdx}
                        onChange={(e) => setCandidateIdx(Number(e.target.value))}
                    >
                        {bundle.candidates.map((c, i) => (
                            <option key={c.original_index} value={i}>
                                #{c.original_index} — {c.affinity != null ? `${c.affinity} kcal/mol` : 'n/a'}
                                {c.default ? ' (best)' : ''}
                            </option>
                        ))}
                    </select>
                    <div className="mt-3 text-xs text-gray-500 font-mono break-all">{candidate.smiles}</div>
                    <div className="mt-3 grid grid-cols-3 gap-2 text-center">
                        <div>
                            <div className="text-lg font-bold text-gray-900">{candidate.affinity ?? '—'}</div>
                            <div className="text-[10px] text-gray-500">kcal/mol</div>
                        </div>
                        <div>
                            <div className="text-lg font-bold text-gray-900">{candidate.qed.toFixed(2)}</div>
                            <div className="text-[10px] text-gray-500">QED</div>
                        </div>
                        <div>
                            <div className="text-lg font-bold text-gray-900">{candidate.sa_score.toFixed(2)}</div>
                            <div className="text-[10px] text-gray-500">SA score</div>
                        </div>
                    </div>
                </div>

                <div className="bg-white p-4 border border-gray-200 shadow-sm rounded-lg">
                    <h3 className="text-sm font-medium text-gray-900 mb-2">Denoising timeline</h3>
                    <div className="flex items-center gap-2 mb-2">
                        <button
                            onClick={() => setPlaying((p) => !p)}
                            className="px-3 py-1 rounded-md text-white text-sm"
                            style={{ backgroundColor: '#065F46' }}
                        >
                            {playing ? 'Pause' : 'Play'}
                        </button>
                        <span className="text-xs text-gray-500 font-mono">Frame {frameIdx} / {nFrames - 1}</span>
                    </div>
                    <input
                        type="range"
                        min={0}
                        max={nFrames - 1}
                        value={frameIdx}
                        onChange={(e) => { setPlaying(false); setFrameIdx(Number(e.target.value)); }}
                        className="w-full"
                    />
                    <div className="flex items-center gap-2 mt-2">
                        <span className="text-xs text-gray-500">Speed</span>
                        <input
                            type="range"
                            min={1}
                            max={10}
                            value={speed}
                            onChange={(e) => setSpeed(Number(e.target.value))}
                            className="w-full"
                        />
                    </div>
                    <p className="text-[11px] text-gray-400 mt-2">
                        Ligand: ball &amp; stick, bonds form as atoms resolve. Receptor pocket
                        (orange) breathes across 6 conformers.
                    </p>
                </div>

                <div className="bg-white p-4 border border-gray-200 shadow-sm rounded-lg">
                    <h3 className="text-sm font-medium text-gray-900 mb-2">Legend</h3>
                    <p className="text-[11px] text-gray-400 mb-2">Hover any atom in the 3D view to see what it is.</p>
                    <div className="grid grid-cols-2 gap-x-3 gap-y-1.5">
                        {ELEMENT_LEGEND.map(({ sym, name, color }) => (
                            <div key={sym} className="flex items-center gap-2 text-xs text-gray-700">
                                <span
                                    className="inline-block w-3 h-3 rounded-full border border-black/10 shrink-0"
                                    style={{ backgroundColor: color }}
                                />
                                <span>{name} <span className="text-gray-400">({sym})</span></span>
                            </div>
                        ))}
                    </div>
                    <div className="mt-3 pt-3 border-t border-gray-100 space-y-1.5">
                        <div className="flex items-center gap-2 text-xs text-gray-700">
                            <span className="inline-block w-3 h-3 rounded-full shrink-0" style={{ backgroundColor: '#22c55e' }} />
                            <span>Ligand carbon (drug candidate)</span>
                        </div>
                        <div className="flex items-center gap-2 text-xs text-gray-700">
                            <span className="inline-block w-3 h-3 rounded-full shrink-0" style={{ backgroundColor: '#f5a623' }} />
                            <span>Pocket-lining residues (binding site)</span>
                        </div>
                        <div className="flex items-center gap-2 text-xs text-gray-700">
                            <span className="inline-block w-3 h-3 rounded-full shrink-0" style={{ backgroundColor: '#5b7fb5' }} />
                            <span>Protein backbone (trypsin)</span>
                        </div>
                    </div>
                </div>
            </div>
        </div>
    );
}
