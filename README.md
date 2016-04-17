# projection-onto-dipole-fields

Implementation of technique described in Liu et al. 2011. A novel background field removal method for MRI using projection onto dipole fields (PDF). NMR Biomedicine.

" For optimal image quality in susceptibility-weighted imaging and accurate quantification of susceptibility, it is necessary to isolate the local field generated by local magnetic sources (such as iron) from the background field that arises from imperfect shimming and variations in magnetic susceptibility of surrounding tissues (including air). A nonparametric background field removal technique based on projection onto dipole fields (PDF). In this PDF technique, the background field inside an ROI is decomposed into a field originating from dipoles outside the ROI using the projection theorem in Hilbert space. This novel PDF background removal technique was validated on a numerical simulation and a phantom experiment and was applied in human brain imaging, demonstrating substantial improvement in background field removal compared with the commonly used high-pass filtering method."

I used this technique to remove the background field in 2D gradient echo (GRE) phase images. 

Input: mask of brain area in 2D, raw phase data with background field.
Output: corrected phase, the field that was removed.

Example image shows the raw phase data, corrected phase, the field that was removed. 
