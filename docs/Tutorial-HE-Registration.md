# Registering an H&E Image

**Prerequisites:** Viewer loaded with a dataset; an H&E image file available (OME-TIFF, standard TIFF, or SVS format)

**Time required:** ~15–30 minutes (depends on landmark placement)

---

## Overview

H&E registration aligns a whole-slide H&E image to the Xenium morphology canvas using a landmark-based similarity transform (rotation, translation, uniform scale). Once registered, the H&E image appears as a layer in napari and you can overlay it against cell segmentations and transcript puncta.

---

## Steps

### 1. Open the H&E tab

In the control panel, open the **Images** group and click the **H&E Registration** tab.

### 2. Load the H&E image

Click **Load H&E Image...** and select your image file. Supported formats include `.tiff`, `.tif`, `.ome.tiff`, and `.svs`.

The image appears as a new layer (`H&E`) in the napari layer list. It may be positioned incorrectly relative to the Xenium image at this stage; that is expected.

<!-- SCREENSHOT: docs/screenshots/tutorial-he-registration-step2.png -->

### 3. Apply flips if needed

If the H&E image appears mirrored or upside down relative to the Xenium morphology image:

- Toggle **Flip vertically** to correct a top-bottom mirror.
- Toggle **Flip horizontally** to correct a left-right mirror.

You can toggle both independently. The canvas updates immediately.

### 4. Run coarse alignment (recommended)

Click **Coarse Align**. The viewer uses tissue outline detection to bring the images into rough alignment automatically. This step is optional but saves time during manual landmark placement.

<!-- SCREENSHOT: docs/screenshots/tutorial-he-registration-step4.png -->

### 5. Place landmark pairs

Landmark pairs link recognisable features in the Xenium image to the same features in the H&E image. You need at least 3 pairs; 5–10 pairs distributed across the tissue give better accuracy.

Good landmark features include:
- Blood vessels (dark lumens in both DAPI and H&E)
- Tissue edges and folds
- Ducts or glands with distinctive shapes

**To place each pair:**

a. Click **Add Xenium Landmark** in the H&E tab. A crosshair cursor appears.
b. Pan and zoom to a recognisable feature in the Xenium morphology image (use the DAPI channel for best contrast). Click once to place the landmark point.
c. Click **Add H&E Landmark**. Navigate to the same feature in the H&E image layer. Click once to place the corresponding point.
d. Repeat steps a–c for each additional feature.

The landmark points appear as coloured dots on the canvas. Matching pairs share the same index number.

<!-- SCREENSHOT: docs/screenshots/tutorial-he-registration-step5.png -->

**Tip:** Switch the active layer in the layer list between `morphology_focus` and `H&E` to compare both images at the same canvas position while placing landmarks.

### 6. Compute the registration

Click **Compute Registration**. The viewer fits a similarity transform to the landmark pairs and applies it to the H&E layer.

The **residuals display** shows per-landmark registration error in micrometres. Typical acceptable values:

| Error | Assessment |
|---|---|
| < 10 µm | Excellent |
| 10–20 µm | Acceptable for tissue-level analysis |
| > 50 µm | Poor; add more landmarks or relocate outliers |

If a specific landmark has a high residual, it may be misplaced. Delete it using the **Remove** button next to its row, re-place it, and recompute.

<!-- SCREENSHOT: docs/screenshots/tutorial-he-registration-step6.png -->

### 7. Inspect the overlay

Use the **H&E opacity** slider in the tab to blend the H&E image against the morphology image. Pan around the tissue at high zoom to check that structures align. Toggle the H&E layer visibility using its eye icon in the layer list to compare with and without the overlay.

### 8. Save the landmarks

Click **Save Landmarks...** and choose a save location. The file is saved as JSON and can be shared with collaborators or used to reproduce the registration on a reinstalled viewer.

<!-- SCREENSHOT: docs/screenshots/tutorial-he-registration-step8.png -->

---

## Notes

- Registration is saved automatically to `sdata_cached.zarr` when you close the viewer. The H&E image path and transform are stored, so the registration is restored on the next launch without repeating these steps.
- To reload a previously saved registration, click **Load Landmarks...** and select the JSON file.
- To start over, click **Clear Landmarks** to remove all landmark points, then repeat from step 5.
- The H&E layer can be used as a reference for drawing annotation polygons; see [Tutorial-Annotations](Tutorial-Annotations).

---

## Next steps

- [Tutorial-ARMS-Overlay](Tutorial-ARMS-Overlay) — register a second image and overlay tile polygons
- [Tutorial-Annotations](Tutorial-Annotations) — draw tissue annotation regions on the registered canvas
- [Tab-HE-Registration](Tab-HE-Registration) — full reference for the H&E Registration tab
