# SMOOTHIE GUI Visual Guide

## Interface Layout

```
┌──────────────────────────────────────────────────────────────────────────┐
│ File   Run   View   Help             [New] [Open] [Save] [Example] [Run] │ ← Menu & Toolbar
├──────────────────────────┬───────────────────────────────────────────────┤
│                          │  ┌─────────────────────────────────────────┐  │
│  INPUT PANEL             │  │  Plot  │  Output Log                    │  │ ← Tab Bar
│                          │  ├─────────────────────────────────────────┤  │
│  ┌────────────────────┐  │  │                                         │  │
│  │ Global   System    │  │  │    [Plot Type: ▼]     [Refresh] [Clear]│  │
│  │ Outgoing Potentials│  │  │                                         │  │
│  ├────────────────────┤  │  │         ┌─────────────────────┐        │  │
│  │                    │  │  │         │                     │        │  │
│  │ Radial Step: 0.05  │  │  │         │                     │        │  │
│  │ Max Radius:  50.0  │  │  │   100 ──┤   Cross Section     │        │  │
│  │                    │  │  │         │                     │        │  │
│  │ Min L:       0     │  │  │    10 ──┤    vs Angle         │        │  │
│  │ Max L:       25    │  │  │         │                     │        │  │
│  │                    │  │  │     1 ──┤                     │        │  │
│  │ Energy:      25.5  │  │  │         │                     │        │  │
│  │ DWBA:        1 ▼   │  │  │   0.1 ──┼─────────────────────┤        │  │
│  │                    │  │  │         0   45   90  135  180         │  │
│  │ Min Angle:   0.0   │  │  │              Angle (deg)               │  │
│  │ Max Angle:   180.0 │  │  │                                         │  │
│  │                    │  │  │                                         │  │
│  │ [Generate Preview] │  │  └─────────────────────────────────────────┘  │
│  │                    │  │                                               │
│  └────────────────────┘  │  OR when "Output Log" tab is selected:       │
│                          │                                               │
│  ← Scrollable input      │  ┌─────────────────────────────────────────┐  │
│     forms with           │  │ [INFO] Starting SMOOTHIE calculation... │  │
│     validation           │  │ [INFO] Reading input file...            │  │
│                          │  │ SMOOTHIE v2.0 - Initializing...         │  │
│                          │  │  Energy = 25.5 MeV                      │  │
│                          │  │  DWBA method = 1                        │  │
│                          │  │  Angular mesh: nx = 34                  │  │
│                          │  │  Radial mesh: nr = 100                  │  │
│                          │  │ Computing cross sections...             │  │
│                          │  │ [SUCCESS] Calculation completed         │  │
│                          │  └─────────────────────────────────────────┘  │
├──────────────────────────┴───────────────────────────────────────────────┤
│ Status: Ready                                           Light Theme [↻] │ ← Status Bar
└──────────────────────────────────────────────────────────────────────────┘
```

## Input Panel - Tab Details

### Tab 1: Global Parameters

```
┌─────────────────────────────────┐
│ ┌──────────────────┐            │
│ │ Global  │ System │            │
│ ├──────────────────┘            │
│ │                               │
│ │ Radial Parameters             │
│ │ ─────────────────             │
│ │ Radial Step Size:  [  0.05 ] │
│ │ Maximum Radius:    [ 50.0  ] │
│ │                               │
│ │ Angular Momentum              │
│ │ ────────────────              │
│ │ Minimum L:         [   0   ] │
│ │ Maximum L:         [  25   ] │
│ │ Minimum Total J:   [   0   ] │
│ │ Maximum Total J:   [  25   ] │
│ │ Minimum lₓ:        [   0   ] │
│ │ Maximum lₓ:        [  12   ] │
│ │                               │
│ │ Energy & Method               │
│ │ ───────────────               │
│ │ Incident Energy:   [ 25.5  ] │
│ │ DWBA Method:       [1: No..▼]│
│ │                               │
│ │ Angular Range                 │
│ │ ─────────────                 │
│ │ Min Angle:         [  0.0  ] │
│ │ Max Angle:         [180.0  ] │
│ │ Increment:         [  1.0  ] │
│ │                               │
│ │ Quadrature                    │
│ │ ──────────                    │
│ │ Gaussian Points:   [  34   ] │
│ │ Radial Points:     [ 100   ] │
│ │                               │
│ │ Output Options                │
│ │ ──────────────                │
│ │ Detailed Output:   [False ▼] │
│ │                               │
│ └───────────────────────────────┘
```

### Tab 2: System Parameters

```
┌─────────────────────────────────┐
│ ┌────────────────────┐          │
│ │ Global │ System    │          │
│ ├────────────────────┘          │
│ │                               │
│ │ ┌─ Projectile (a) ─────────┐ │
│ │ │ Name:         [  d    ]  │ │
│ │ │ Mass (amu):   [ 2.0   ]  │ │
│ │ │ Charge:       [ 1.0   ]  │ │
│ │ │ Spin:         [ 0.0   ]  │ │
│ │ └──────────────────────────┘ │
│ │                               │
│ │ ┌─ Target (A) ──────────────┐ │
│ │ │ Name:         [ 93Nb  ]  │ │
│ │ │ Mass (amu):   [ 93.0  ]  │ │
│ │ │ Charge:       [ 41.0  ]  │ │
│ │ │ Spin:         [ 0.0   ]  │ │
│ │ └──────────────────────────┘ │
│ │                               │
│ │ ┌─ Detected Particle (b) ──┐ │
│ │ │ Name:         [  p    ]  │ │
│ │ │ Mass (amu):   [1.0078 ]  │ │
│ │ │ Charge:       [ 1.0   ]  │ │
│ │ │ Spin:         [ 0.0   ]  │ │
│ │ └──────────────────────────┘ │
│ │                               │
│ │ ┌─ Undetected Particle (x) ┐ │
│ │ │ Name:         [  n    ]  │ │
│ │ │ Mass (amu):   [1.0087 ]  │ │
│ │ │ Charge:       [ 0.0   ]  │ │
│ │ │ Spin:         [ 0.0   ]  │ │
│ │ └──────────────────────────┘ │
│ │                               │
│ │ ┌─ System Properties ───────┐ │
│ │ │ Spin Coupling:  [ 0.0 ]  │ │
│ │ │ Orbital L:      [  0  ]  │ │
│ │ │ Nodes:          [  1  ]  │ │
│ │ │ Binding (MeV):  [2.224]  │ │
│ │ └──────────────────────────┘ │
│ └───────────────────────────────┘
```

### Tab 3: Outgoing Parameters

```
┌─────────────────────────────────┐
│ ┌───────────────────────┐       │
│ │ Global │ System │     │       │
│ │ Outgoing │ Potentials │       │
│ ├───────────────────────┘       │
│ │                               │
│ │ Energy Range                  │
│ │ ────────────                  │
│ │                               │
│ │ Minimum Energy: [  2.0  ] MeV│
│ │                               │
│ │ Maximum Energy: [ 30.0  ] MeV│
│ │                               │
│ │ Energy Step:    [  1.0  ] MeV│
│ │                               │
│ │                               │
│ └───────────────────────────────┘
```

### Tab 4: Potentials

```
┌─────────────────────────────────┐
│ ┌──────────────────────────┐    │
│ │ Global │ System │ Outgoing│   │
│ │ Potentials                │    │
│ ├──────────────────────────┘    │
│ │                               │
│ │ ┌────────────────────────┐   │
│ │ │Pot 1│Pot 2│Pot 3│Pot 4│   │
│ │ │ (a) │ (b) │ (x) │ (p) │   │
│ │ ├────────────────────────┘   │
│ │ │                             │
│ │ │ Potential 1: a (proj-tgt)  │
│ │ │ ─────────────────────────  │
│ │ │                             │
│ │ │ Model:    [Woods-Saxon ▼]  │
│ │ │                             │
│ │ │ Mass 1:   [   0.0  ]  amu  │
│ │ │ Mass 2:   [  93.0  ]  amu  │
│ │ │ Coulomb R:[   1.3  ]  fm   │
│ │ │                             │
│ │ │ Real Volume                 │
│ │ │ ───────────                 │
│ │ │ Depth:    [  77.3  ]  MeV  │
│ │ │ Radius:   [   1.15 ]  fm   │
│ │ │ Diffuse:  [   0.77 ]  fm   │
│ │ │                             │
│ │ │ Imaginary Volume            │
│ │ │ ────────────────            │
│ │ │ Depth:    [   6.1  ]  MeV  │
│ │ │ Radius:   [   1.33 ]  fm   │
│ │ │ Diffuse:  [   0.47 ]  fm   │
│ │ │                             │
│ │ │ Surface Absorption          │
│ │ │ ──────────────────          │
│ │ │ Depth:    [   8.4  ]  MeV  │
│ │ │ Radius:   [   1.37 ]  fm   │
│ │ │ Diffuse:  [   0.77 ]  fm   │
│ │ │                             │
│ │ └─────────────────────────────┘
│ └───────────────────────────────┘
```

## Plot Widget - Different Views

### Cross Section vs Angle (Default)

```
┌──────────────────────────────────────┐
│ Plot Type: [Cross Section vs Angle▼]│
│ [Refresh] [Clear]                    │
├──────────────────────────────────────┤
│                                      │
│  d$\sigma$/d$\Omega$ (mb/sr)         │
│  │                                   │
│  │  ●                                │
│10²│   ●                               │
│  │    ●                              │
│  │     ●●                            │
│10¹│       ●●                          │
│  │         ●●●                       │
│  │            ●●●                    │
│10⁰│               ●●●●                │
│  │                   ●●●●●           │
│  │                        ●●●●●●     │
│10⁻¹│                            ●●●●●│
│  └────────────────────────────────── │
│    0   30   60   90  120  150  180  │
│           Angle (degrees)            │
│                                      │
│  [🔍] [⌂] [◀] [▶] [⊕] [⊖] [💾]     │ ← Matplotlib toolbar
└──────────────────────────────────────┘
```

### Angular Distribution (Polar)

```
┌──────────────────────────────────────┐
│ Plot Type: [Angular Distribution  ▼]│
│ [Refresh] [Clear]                    │
├──────────────────────────────────────┤
│                                      │
│            0°                        │
│            ●                         │
│          ●   ●                       │
│        ●       ●                     │
│      ●           ●                   │
│ 270°●             ●  90°             │
│      ●           ●                   │
│        ●       ●                     │
│          ●   ●                       │
│            ●                         │
│           180°                       │
│                                      │
│  Polar plot showing angular          │
│  distribution of cross section       │
│                                      │
│  [🔍] [⌂] [◀] [▶] [⊕] [⊖] [💾]     │
└──────────────────────────────────────┘
```

## Color Scheme

### Light Theme
- Background: #ffffff (white)
- Text: #333333 (dark gray)
- Accent: #007AFF (blue)
- Borders: #e0e0e0 (light gray)
- Buttons: #007AFF → #0051D5 (hover)

### Dark Theme
- Background: #2d2d2d (dark gray)
- Text: #e0e0e0 (light gray)
- Accent: #0A84FF (bright blue)
- Borders: #404040 (medium gray)
- Buttons: #0A84FF → #0070E0 (hover)

## Status Indicators

### Log Colors

```
[INFO]    Blue:   #0066CC  →  Informational messages
[WARNING] Orange: #FF8800  →  Warnings, non-critical issues
[ERROR]   Red:    #CC0000  →  Errors, calculation failed
[SUCCESS] Green:  #009900  →  Successful completion
Output    Black:  #333333  →  SMOOTHIE stdout/stderr
```

### Progress States

```
Ready           →  Idle, waiting for user action
Running...      →  Calculation in progress
Completed ✓     →  Success
Failed ✗        →  Error occurred
Stopped         →  User interrupted
```

## Responsive Design

### Window Sizes

**Minimum**: 1200 × 700 px
**Recommended**: 1600 × 900 px
**Optimal**: 1920 × 1080 px

### Splitter Behavior

```
Small (1200px):   40% Input | 60% Output
Medium (1600px):  35% Input | 65% Output
Large (1920px):   30% Input | 70% Output
```

User can drag splitter to adjust ratios.

## Keyboard Navigation

### Global Shortcuts

```
Ctrl/Cmd + N     →  New file
Ctrl/Cmd + O     →  Open file
Ctrl/Cmd + S     →  Save file
Ctrl/Cmd + R     →  Run calculation
Ctrl/Cmd + .     →  Stop calculation
Ctrl/Cmd + Q     →  Quit application
Ctrl/Cmd + F     →  Find in log (when focused)
Tab              →  Navigate between fields
Enter            →  Confirm/next field
```

### In-Form Navigation

```
Tab              →  Next field
Shift + Tab      →  Previous field
↑ ↓              →  Increment/decrement spinboxes
Space            →  Toggle checkboxes
Enter            →  Activate buttons
```

## Example Workflow Diagram

```
┌─────────────┐
│   Launch    │
│     GUI     │
└──────┬──────┘
       │
       ▼
┌─────────────┐     ┌─────────────┐
│ Load Example│────▶│Review Params│
└─────────────┘     └──────┬──────┘
       │                   │
       │                   ▼
       │            ┌─────────────┐
       │            │Modify Values│
       │            └──────┬──────┘
       │                   │
       └───────────────────┘
                   │
                   ▼
            ┌─────────────┐
            │  Save Input │ (optional)
            └──────┬──────┘
                   │
                   ▼
            ┌─────────────┐
            │Run Smoothie │
            └──────┬──────┘
                   │
                   ▼
            ┌─────────────┐
            │Monitor Log  │
            └──────┬──────┘
                   │
                   ▼
            ┌─────────────┐
            │View Plots   │
            └──────┬──────┘
                   │
                   ▼
            ┌─────────────┐
            │Export Results│
            └─────────────┘
```

---

This visual guide shows the complete interface layout and interaction patterns of the SMOOTHIE GUI. For detailed usage instructions, see `README.md` and `GUI_QUICKSTART.md`.
