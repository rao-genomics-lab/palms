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

### Automatic fine alignment

| Control | Description |
|---|---|
| Fine Align (nuclei) | Refines the coarse transform by matching every nucleus in the H&E to Xenium's own nuclear masks. Reports how many nuclei matched, the median residual, and how far it moved the image. Disabled until a coarse (or landmark) transform exists, and on a dataset with no `nucleus_labels`. |

This is what the manual landmark step is doing, done over the whole section
instead of the handful of points you can click: haematoxylin is deconvolved out
of the H&E, its nuclei are found as sub-pixel peaks, and that point set is fitted
onto the centroids of `labels/nucleus_labels` — Xenium's segmentation of the DAPI
channel. On the two reference datasets it takes the coarse transform's 15–17 µm
down to **0.7–0.8 µm**, with nothing placed by hand. That is better than the
landmark fit reaches, and better than 10x's own shipped matrix: on nuclei held
out of the fit entirely, the automatic transform places them closer to the
nuclear masks than either.

It needs a starting transform, which is what Coarse Align is for, so run that
first. It corrects the seed in two passes — a search over scale and rotation,
then the per-nucleus fit — which between them recover from a coarse alignment
60 µm out. Confidence is an **enrichment over chance**: how many more H&E nuclei
land within a micron of a nuclear mask than a random point set of the same
density would. A real fit scores 9–22×; a transform a few degrees wrong, or an
H&E of a different section, scores about 1×.

If it reports **LOW CONFIDENCE**, the usual cause is the starting transform, not
the H&E and not the section. The message shows how much of the move came from the
scale/rotation search; a large figure there means Coarse Align's scale was
materially out, which happens most on big sections where its tissue-outline scale
prior is weakest. Re-running Coarse Align, or placing three landmarks and
pressing "Compute Registration" to give the nuclei fit a better start, both work.

It takes one to three minutes on a whole section, most of it finding nuclei in
the full-resolution H&E. That resolution is the point: the same detector run one
pyramid level down is three times less accurate, and no number of extra
detections makes up for it, because a coarser pixel biases each peak rather than
adding noise that averages away.

You can still place landmarks afterwards — "Compute Registration" overwrites the
automatic fit — but on a dataset where the nuclei fit reports high confidence
there is usually nothing left to correct by hand.

### Landmark registration

| Control | Description |
|---|---|
| Add Xenium Landmark | Activates the Xenium landmark layer in "add point" mode. Click a recognisable feature in the Xenium image to place a landmark. |
| Add H&E Landmark | Activates the H&E landmark layer in "add point" mode. Click the corresponding feature in the H&E image. |
| Clear All | Removes all landmarks and clears the fine registration affine, including an automatic one. |
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
4. Click "Fine Align (nuclei)". This is normally the last step — it refines the
   coarse transform to well under a micron by matching the H&E's nuclei to the
   nuclear masks, and takes one to three minutes. Check the reported enrichment;
   if it says LOW CONFIDENCE, fall back to landmarks below.
5. Only if the automatic fit was refused or looks wrong — place at least 3
   landmark pairs:
   1. Click "Add Xenium Landmark", then click a recognisable feature (for example, a vessel or duct) in the Xenium morphology image.
   2. Click "Add H&E Landmark", then click the same feature in the H&E image.
   3. Repeat for at least 3 features spread across the tissue.
6. Click "Compute Registration". Inspect the residuals — values under 10–20 µm
   are typically acceptable. This **replaces** an automatic fit.
7. (Optional) Click "Save Landmarks..." to save the result for sharing or future reuse.

## Notes

- "Save Landmarks..." writes the H&E points in the orientation they were fitted
  in, together with the flip settings, so the file is self-consistent: its affine
  maps its own H&E points onto its own Xenium points. `scripts/compare_he_registration.py`
  scores such a file against 10x's shipped alignment matrix, and
  `scripts/score_coarse_align.py` scores Coarse Align against the same reference.
- **A saved `landmarks.json` is worth keeping.** The fine transform has one slot,
  written by whichever method ran last, so an automatic fit replaces a landmark
  fit in the session. "Save Landmarks..." writes a record that survives, and
  `scripts/score_nuclei_align.py --landmarks` will score against it.
- `scripts/score_nuclei_align.py` scores the whole automatic chain the same way,
  and adds the check that agreement with another estimate cannot give: it fits on
  one half of the section and scores every candidate transform on the nuclei of
  the other half, which that fit has never seen. All three estimates — 10x's,
  the landmark fit and the automatic one — agree with each other only to about a
  micron, while the automatic fit reproduces itself to 0.02 µm, so at this
  accuracy the reference is no longer a fine enough ruler and the held-out
  residual is what to read.
- Registration is persisted automatically to `sdata_cached.zarr` and restored on the next launch; you do not need to redo registration between sessions. The transform lives in `viewer_session/he/`, alongside the image path, its dimensions and the flip settings.
- The H&E image itself is also stored in the zarr cache, so you do not need to re-load the original file on subsequent sessions.
- The flips, the coarse alignment, the automatic nuclei fit and the landmark registration are each recorded into the analysis provenance, with the landmark coordinates inlined, so the exported notebook reproduces the transform without needing the landmark files.
