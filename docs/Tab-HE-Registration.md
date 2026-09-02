# H&E Registration

Load an H&E (or other brightfield) TIFF image, align it to the Xenium coordinate system using optional coarse tissue-outline alignment followed by manual landmark-based registration, and persist the affine transform. This tab is in the "Images" control panel group.

![He Registration](screenshots/tab-he-registration.png)

## Controls

### Loading and orientation

| Control | Description |
|---|---|
| Load H&E Image... | Opens a file dialog accepting `.ome.tif`, `.tif`, `.tiff`, `.svs` files. Loads the image as a multi-scale pyramid and adds it to the napari canvas. |
| Flip vertically | Flips the H&E image vertically by applying a flip affine. |
| Flip horizontally | Flips the H&E image horizontally. |
| Opacity | Slider (0–100, default 70). Controls H&E layer transparency. Disabled until an image is loaded. |

### Coarse alignment

| Control | Description |
|---|---|
| Coarse Align | Automatically overlays the H&E on the Xenium morphology image, by searching rotation, scale and reflection for the best match between blurred nuclear density in the two modalities. Reports the scale, rotation, whether the H&E is mirrored, and a match score; a **mirrored H&E ticks "Flip horizontally" itself**. Disabled until an image is loaded. |

Coarse Align takes a few seconds — it scores a few thousand orientations rather
than guessing one — and it says how much to believe the result. A low match
score, or one that barely beats the next distinct orientation, is reported as
**LOW CONFIDENCE**: the transform is still applied, but check the overlay before
placing landmarks. On the two reference datasets it lands within 0.7% of the
correct scale and 0.04° of the correct rotation, which is 15–17 µm across the
whole slide.

Two things affect how well it does. If the H&E declares its own pixel size
(most OME-TIFFs do), that is used as the scale prior and is worth about 0.06%;
otherwise the scale is estimated from the two tissue outlines. And if tissue
fills both images — which is what a **crop export** looks like — the outlines
carry no information, so the match is made on internal structure instead. That
case works, but scores lower than a whole section does.

### Landmark registration

| Control | Description |
|---|---|
| Add Xenium Landmark | Activates the Xenium landmark layer in "add point" mode. Click a recognisable feature in the Xenium image to place a landmark. |
| Add H&E Landmark | Activates the H&E landmark layer in "add point" mode. Click the corresponding feature in the H&E image. |
| Clear All | Removes all landmarks and clears the fine registration affine. |
| Compute Registration | Computes a similarity affine from the placed landmark pairs. Enabled when 3 or more pairs are present. Displays per-landmark residuals (pixels and µm), mean and max residuals, and the computed scale factor. |
| Residuals (read-only) | Text area showing registration quality metrics after the last computation. |
| Save Landmarks... | Saves landmarks and the computed affine to a JSON file. Enabled when at least one landmark is present. |
| Load Landmarks... | Restores landmarks and affine from a previously saved JSON file. |

The combined affine applied to the H&E image is composed as the flip affine
followed by the landmark affine, or by the coarse affine when no landmarks have
been placed. Both are fitted in the *flipped* frame, so toggling a flip
afterwards clears a coarse alignment — it was fitted for the other orientation
— and you should run Coarse Align again.

## Workflow

1. Click "Load H&E Image..." and select your file.
2. Click "Coarse Align". It finds the rotation, scale and reflection itself, so
   you do not need to set the flips by hand first — it ticks "Flip horizontally"
   for you if the H&E turns out to be mirrored. Check the reported match score.
3. If the overlay is wrong, or Coarse Align reports low confidence, set "Flip
   vertically" / "Flip horizontally" yourself and run it again, or go straight
   to landmarks.
4. Place at least 3 landmark pairs:
   1. Click "Add Xenium Landmark", then click a recognisable feature (for example, a vessel or duct) in the Xenium morphology image.
   2. Click "Add H&E Landmark", then click the same feature in the H&E image.
   3. Repeat for at least 3 features spread across the tissue.
5. Click "Compute Registration". Inspect the residuals — values under 10–20 µm are typically acceptable.
6. (Optional) Click "Save Landmarks..." to save the result for sharing or future reuse.

## Notes

- "Save Landmarks..." writes the H&E points in the orientation they were fitted
  in, together with the flip settings, so the file is self-consistent: its affine
  maps its own H&E points onto its own Xenium points. `scripts/compare_he_registration.py`
  scores such a file against 10x's shipped alignment matrix, and
  `scripts/score_coarse_align.py` scores Coarse Align against the same reference.
- Registration is persisted automatically to `sdata_cached.zarr` and restored on the next launch; you do not need to redo registration between sessions. The transform lives in `viewer_session/he/`, alongside the image path, its dimensions and the flip settings.
- The H&E image itself is also stored in the zarr cache, so you do not need to re-load the original file on subsequent sessions.
- The flips, the coarse alignment and the landmark registration are each recorded into the analysis provenance, with the landmark coordinates inlined, so the exported notebook reproduces the transform without needing the landmark files.
