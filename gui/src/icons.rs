//! Remix Icon constants for use in the GUI
//! Full icon list: https://remixicon.com/
//!
//! Usage: ui.button(format!("{} Save", icons::SAVE));

// === File Operations ===
pub const FOLDER_OPEN: &str = "\u{ed70}"; // ri-folder-open-line
pub const FOLDER_ADD: &str = "\u{ed5a}"; // ri-folder-add-line
pub const FILE: &str = "\u{ecc3}"; // ri-file-2-line
pub const FILE_ADD: &str = "\u{ecc9}"; // ri-file-add-line
pub const SAVE: &str = "\u{f0b3}"; // ri-save-line
pub const DOWNLOAD: &str = "\u{ec5a}"; // ri-download-line
pub const UPLOAD: &str = "\u{f25c}"; // ri-upload-line
pub const IMPORT: &str = "\u{ee54}"; // ri-import-line (actually login)
pub const EXPORT: &str = "\u{ee96}"; // ri-share-line

// === Playback / Simulation ===
pub const PLAY: &str = "\u{f00b}"; // ri-play-line
pub const PLAY_CIRCLE: &str = "\u{f009}"; // ri-play-circle-line
pub const PAUSE: &str = "\u{efd8}"; // ri-pause-line
pub const STOP: &str = "\u{f1a1}"; // ri-stop-line
pub const STOP_CIRCLE: &str = "\u{f19f}"; // ri-stop-circle-line
pub const LOADER: &str = "\u{eeca}"; // ri-loader-line
pub const REFRESH: &str = "\u{f064}"; // ri-refresh-line

// === View / Display ===
pub const EYE: &str = "\u{ecb5}"; // ri-eye-line
pub const EYE_OFF: &str = "\u{ecb7}"; // ri-eye-off-line
pub const FULLSCREEN: &str = "\u{ed9c}"; // ri-fullscreen-line
pub const FULLSCREEN_EXIT: &str = "\u{ed9a}"; // ri-fullscreen-exit-line
pub const ZOOM_IN: &str = "\u{f2db}"; // ri-zoom-in-line
pub const ZOOM_OUT: &str = "\u{f2dd}"; // ri-zoom-out-line
pub const FOCUS: &str = "\u{ed4e}"; // ri-focus-line (fit/center view)
pub const GRID: &str = "\u{eddf}"; // ri-grid-line
pub const LAYOUT: &str = "\u{ee95}"; // ri-layout-line

// === 3D View Controls ===
pub const ROTATE_3D: &str = "\u{f08a}"; // ri-rotate-lock-line (3D rotation)
pub const FLIP_H: &str = "\u{ed24}"; // ri-flip-horizontal-line
pub const FLIP_V: &str = "\u{ed26}"; // ri-flip-vertical-line
pub const CONTRAST: &str = "\u{ebd0}"; // ri-contrast-line
pub const ASPECT_RATIO: &str = "\u{ea7e}"; // ri-aspect-ratio-line

// === View Directions (using arrow boxes) ===
pub const VIEW_FRONT: &str = "\u{ea5f}"; // ri-arrow-left-s-line (pointing at viewer)
pub const VIEW_BACK: &str = "\u{ea6d}"; // ri-arrow-right-s-line
pub const VIEW_LEFT: &str = "\u{ea77}"; // ri-arrow-up-s-line
pub const VIEW_RIGHT: &str = "\u{ea4b}"; // ri-arrow-down-s-line
pub const VIEW_TOP: &str = "\u{eab8}"; // ri-arrow-up-circle-line
pub const VIEW_BOTTOM: &str = "\u{ea9a}"; // ri-arrow-down-circle-line
pub const VIEW_ISO: &str = "\u{eae4}"; // ri-box-3-line (isometric)

// === Camera Controls ===
pub const CAMERA_RESET: &str = "\u{f064}"; // ri-refresh-line (reset camera)
pub const CAMERA_SWITCH: &str = "\u{eb37}"; // ri-camera-switch-line
pub const CROSSHAIR: &str = "\u{ed48}"; // ri-focus-2-line (center point)

// === Clipping / Sectioning ===
pub const SCISSORS: &str = "\u{f0c1}"; // ri-scissors-line
pub const CUT: &str = "\u{f0c1}"; // ri-scissors-line (alias)
pub const SLICE: &str = "\u{ec5e}"; // ri-divide-line
pub const SPLIT: &str = "\u{f17b}"; // ri-split-cells-horizontal

// === Display Modes ===
pub const WIREFRAME: &str = "\u{eddf}"; // ri-grid-line (wireframe)
pub const SHADED: &str = "\u{ee9e}"; // ri-lightbulb-line (shaded)
pub const EDGES: &str = "\u{f0c3}"; // ri-scan-line (edge display)
pub const FACES: &str = "\u{f3cd}"; // ri-hexagon-fill (polygon face)
pub const POINTS: &str = "\u{eb7c}"; // ri-checkbox-blank-circle-fill (node dot)

// === Annotations ===
pub const STICKY_NOTE: &str = "\u{f198}"; // ri-sticky-note-line
pub const MARKUP: &str = "\u{ef18}"; // ri-markup-line
pub const TEXT: &str = "\u{f1fa}"; // ri-text
pub const PENCIL: &str = "\u{efc6}"; // ri-pencil-line
pub const HIGHLIGHT: &str = "\u{ef12}"; // ri-mark-pen-line

// === Measurement ===
pub const RULER_2: &str = "\u{f0a1}"; // ri-ruler-2-line
pub const DISTANCE: &str = "\u{f04d}"; // ri-radar-line (measure)

// === Edit / Transform ===
pub const EDIT: &str = "\u{ec86}"; // ri-edit-line
pub const DELETE: &str = "\u{ec2a}"; // ri-delete-bin-line
pub const ADD: &str = "\u{ea13}"; // ri-add-line
pub const ADD_CIRCLE: &str = "\u{ea11}"; // ri-add-circle-line
pub const SUBTRACT: &str = "\u{f1a7}"; // ri-subtract-line
pub const DRAG_MOVE: &str = "\u{ec62}"; // ri-drag-move-line
pub const CURSOR: &str = "\u{ec0a}"; // ri-cursor-line

// === Settings / Config ===
pub const SETTINGS: &str = "\u{f0ee}"; // ri-settings-line
pub const SETTINGS_4: &str = "\u{f0e8}"; // ri-settings-4-line (gear)
pub const EQUALIZER: &str = "\u{ec9d}"; // ri-equalizer-line (sliders)
pub const TOOLS: &str = "\u{f222}"; // ri-tools-line
pub const HAMMER: &str = "\u{edef}"; // ri-hammer-line

// === Shapes / Geometry (for FEA) ===
pub const BOX_3: &str = "\u{eae4}"; // ri-box-3-line (3D box)
pub const SHAPE: &str = "\u{f0df}"; // ri-shape-line
pub const ARTBOARD: &str = "\u{ea7c}"; // ri-artboard-line
pub const RULER: &str = "\u{f0a1}"; // ri-ruler-2-line
pub const COMPASS: &str = "\u{ebc4}"; // ri-compass-line
pub const NODE_TREE: &str = "\u{ef90}"; // ri-node-tree

// === Info / Status ===
pub const INFO: &str = "\u{ee59}"; // ri-information-line
pub const QUESTION: &str = "\u{f045}"; // ri-question-line
pub const WARNING: &str = "\u{eca1}"; // ri-error-warning-line
pub const CHECK: &str = "\u{eb7b}"; // ri-check-line
pub const CHECK_CIRCLE: &str = "\u{eb81}"; // ri-checkbox-circle-line
pub const CLOSE: &str = "\u{eb99}"; // ri-close-line
pub const CLOSE_CIRCLE: &str = "\u{eb97}"; // ri-close-circle-line

// === Navigation / Menu ===
pub const MENU: &str = "\u{ef3e}"; // ri-menu-line
pub const HOME: &str = "\u{ee2b}"; // ri-home-line
pub const ARROW_LEFT: &str = "\u{ea60}"; // ri-arrow-left-line
pub const ARROW_RIGHT: &str = "\u{ea6e}"; // ri-arrow-right-line
pub const ARROW_UP: &str = "\u{ea78}"; // ri-arrow-up-line
pub const ARROW_DOWN: &str = "\u{ea4c}"; // ri-arrow-down-line
pub const EXTERNAL_LINK: &str = "\u{ecaf}"; // ri-external-link-line

// === Screenshot / Image ===
pub const SCREENSHOT: &str = "\u{f0c7}"; // ri-screenshot-line
pub const CAMERA: &str = "\u{eb31}"; // ri-camera-line
pub const IMAGE: &str = "\u{ee4b}"; // ri-image-line
pub const GALLERY: &str = "\u{eda5}"; // ri-gallery-line

// === Misc ===
pub const CLIPBOARD: &str = "\u{eb91}"; // ri-clipboard-line
pub const TIME: &str = "\u{f20a}"; // ri-time-line
pub const HISTORY: &str = "\u{ee14}"; // ri-history-line
