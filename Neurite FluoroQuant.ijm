// Neurite FluoroQuant - Fast Compute — same measurements as Neurite FluoroQuant,
// with ImageJ batch mode hiding non-interactive windows and one-pass intensity scans.
// Install: copy to Fiji.app/plugins/ then Help > Refresh Menus.
// Menu: Plugins > Neurite FluoroQuant - Fast Compute
// Does not replace Neurite FluoroQuant.ijm.
//
MACRO_VERSION = "2026-11-05-fast-compute";
//
// Traces a neurite on the ROI channel, straightens both channels along that path,
// builds a binary mask from the ROI channel, and writes optional intensity, CV,
// neuropil, and patch measurements. Outputs are CSV/PNG files under Analysis/,
// ROI/, and Parameters/ in a numbered run folder.
//
// Key features:
//  * User-assigned ROI and analysis channels; mask from the ROI channel along the trace.
//  * 2D, Z-projection, or per-plane 3D tracing, with optional background subtraction.
//  * Optional in-mask intensity, CV (fluorescence and width), neuropil, and patch detection.
//  * Optional adjacent-particle mode for signal beside the neurite (not mask-gated),
//    with midline-distance filtering and an optional second analysis channel.
//  * Interactive tuning for mask threshold and patch detection; batch reuse of saved settings.

requires("1.52o15");

// Global work buffers for region detection (ImageJ limits user functions to 20 arguments).
var peakRegAreaWork, peakRegMeanWork;
var peakRegBXWork, peakRegBYWork, peakRegWWork, peakRegHWork, peakRegFeretWork, peakRegMinFeretWork, peakRegCircWork;
var peakRegionKeptRoiIdx;
// Region separation tuning (passed via globals; runPeakDetection is already at the 20-argument limit).
var regionSepMode, regionPromPct, regionDogUm;
// Region stage counts for tuning dialog (threshold-only vs post-watershed).
var regionCountThr, regionCountThrFrag, regionCountWs, regionCountFrag, regionResidualSpan;
// Particles that passed min area, before the adjacent midline-distance filter.
var regionCountPreDist;
// Sub-min-area speckle outlines for diagnostic preview when final count is zero.
var regionDiagN, regionDiagBXWork, regionDiagBYWork, regionDiagWWork, regionDiagHWork;
// Saved ROI set (true pixel outlines) of the diagnostic speckles, so the orange preview
// can trace actual shapes instead of bounding boxes. Empty path = not available.
var regionDiagRoiPath = "";
// Saved ROI set of kept regions for final overlay PNG export (true pixel outlines).
var regionExportRoiPath = "";
// 3D per-plane path tracing state: one ROI per segment in the ROI Manager,
// segZ[i] = slice index for segment i, pathTotalPx = summed in-plane arc length.
var segZ = newArray(0);
var segStartX = newArray(0);
var segStartY = newArray(0);
var segEndX = newArray(0);
var segEndY = newArray(0);
var pathSegCount = 0;
var pathTotalPx = 0;
var projZFirst = 0;
var projZLast = 0;
var zProjMethod = "Max Intensity";
var compensateZIntensity = true;
var batchMode = 0;
var batchProcessedCount = 0;
var roiIdx, ch2Idx, fovUm, straightenUm, tissueType;
var compensateZIntensityPref, is3DProject, is3DPerPlane, is3D;
var maskIdx, rbUm, maskThr, maskDilUm;
var includeCh2IntSummary, includeCh2Profile, includeNP;
var includeCvSummary, includeCv;
var includeRoiSummary, includeRoiProfile, includeMaskWidth, includeAnalysisParams;
var roiCvGlobal, maskWidthCv;
var scanSumIn, scanCntIn, scanSumOut, scanCntOut;
var includeCh1, wcvUm, doDetrend, detrendUm;
var maskTitle, analysisCorrWin;
var idAnalysisCorr, idAnalysis2Corr, idMaskCore;
var hiddenComputeOn;
var analysisExportTag, analysisOverlayWin, peakPreviewPrefix, peakProfilePlotWin, analysisPlainWin;
var analysisOutDir, analysisWin, overlayPeaksPath;
var swPk, shPk, umPerPx, roiLenUm, base, outDir, autoIncrementOutput;
var widthsUm, widthsPx, maskWidthMeanUm, totalMaskArea, peakDensity;
var tolPct, minSepUm, markerUm, doRefine;
var regionMinAreaUm, regionMaxAreaUm, regionSmoothUm;
var peaksBright, peakPreviewZoomPct;
var def_tolPct, def_minSepUm, def_markerUm, def_doRefine;
var def_regionMinAreaUm, def_regionMaxAreaUm, def_regionSmoothUm, def_regionPromPct, def_regionDogUm;
var def_adj2_tolPct, def_adj2_minSepUm, def_adj2_regionMinAreaUm, def_adj2_regionMaxAreaUm;
var def_adj2_regionSmoothUm, def_adj2_regionPromPct, def_adj2_regionDogUm;
var def_adjMaxDistMidUm, def_adjSecondAnalysisIdx;
var swArc, shArc, profileArc, yAtMax, arcValid;
var profileArcExport, arcValidExport, swArcExport;
var arcProfMin, arcProfMax, arcProfValidCols;
var nPeaks, nPeaksArcFM, nPeaksArcSep;
var doPeakDetection, isRegionMethodUsed, detectionType, peakMethod;
var xs, ys, vals;
var regAreaExport, regMeanExport;
var regBXExport, regBYExport, regWExport, regHExport, regFeretExport, regMinFeretExport, regCircExport;
var appendMasterSummary, masterSummaryPath;
// Adjacent-particle mode: detect fluorescence outside the mask (e.g. synaptic contacts onto a
// dendrite). The mask is still built, but only as the reference edge for per-particle distances.
var adjacentMode = 0;
// Optional second synaptic/analysis channel in adjacent mode only (0 = none). Within-mask mode ignores it.
var adjSecondAnalysisIdx = 0;
var adjDualMode = 0;
// Analysis channel index for the active adjacent pass (tuning titles / ROI zip names).
var masterAnalysisChIdx = 0;
// 1 while running the second adjacent detection pass (use unified.adj2.* prefs).
var adjPassUseAdj2Prefs = 0;
// Second analysis channel windows (adjacent dual only).
var analysis2Win, analysis2StraightWin, analysis2CorrWin, analysis2PlainWin, analysis2OverlayWin;
var analysis2ExportTag;
var analysis2MaskPreviewWin = "";
var analysis2BgsubPreviewWin = "";
// ROI grayscale kept for dual-adjacent merged RGB overlay (includeCh1 is off in adjacent mode).
var adjMergeRoiWin = "";
var adjRoiZipPass1 = "";
var adjRoiZipPass2 = "";
var adjMergedOverlayWin = "Adjacent_merged_overlay";
// Stashed adjacent master-CSV metrics (one row per run; pass 1 = primary, pass 2 = second channel).
var adjMasterReady = 0;
var adjMasterChIdx1 = 0;
var adjMasterNPeaks1 = 0;
var adjMasterPeakDensity1 = 0;
var adjMasterMeanArea1 = 0;
var adjMasterMeanInt1 = 0;
var adjMasterMeanFeret1 = 0;
var adjMasterMeanMinFeret1 = 0;
var adjMasterMeanCirc1 = 0;
var adjMasterMeanEquiv1 = 0;
var adjMasterLinear1 = 0;
var adjMasterDistMask1 = 0;
var adjMasterDistMid1 = 0;
var adjMasterChIdx2 = 0;
var adjMasterNPeaks2 = 0;
var adjMasterPeakDensity2 = 0;
var adjMasterMeanArea2 = 0;
var adjMasterMeanInt2 = 0;
var adjMasterMeanFeret2 = 0;
var adjMasterMeanMinFeret2 = 0;
var adjMasterMeanCirc2 = 0;
var adjMasterMeanEquiv2 = 0;
var adjMasterLinear2 = 0;
var adjMasterDistMask2 = 0;
var adjMasterDistMid2 = 0;
// Stashed adjacent detection tuning for parameters CSV (Analysis Ch1 / Ch2).
var adjParamTol1 = 0;
var adjParamMinSep1 = 0;
var adjParamMinArea1 = 0;
var adjParamMaxArea1 = 0;
var adjParamSmooth1 = 0;
var adjParamProm1 = 0;
var adjParamDog1 = 0;
var adjParamTol2 = 0;
var adjParamMinSep2 = 0;
var adjParamMinArea2 = 0;
var adjParamMaxArea2 = 0;
var adjParamSmooth2 = 0;
var adjParamProm2 = 0;
var adjParamDog2 = 0;
// Set after adjacent per-channel exports so the shared export block does not duplicate them.
var adjacentExportsDone = 0;
// Particles farther than this from the straightened midline are rejected. Distances to the mask are
// still measured and reported, but gating on them cannot separate the neurite from a structure
// running alongside it, so the midline is the only criterion.
var adjMaxDistMidUm = 1.5;
// Mask component that sits on the straightened midline, i.e. the traced neurite itself. Off-neurite
// thresholded fragments would otherwise sit within the gating distance of unrelated particles.
var adjCoreMaskTitle = "";
// Pixel size of the straightened band, so the preview can draw the acceptance band in microns.
var adjUmPerPx = 0;
// Per-particle distances (um) and the originating ROI Manager index, parallel to peakReg*Work.
var peakRegDistMaskWork, peakRegDistMaskCentWork, peakRegDistMidWork, peakRegSrcIdxWork;
var regDistMaskExport, regDistMaskCentExport, regDistMidExport;
var adjLinearDensity = 0;
// Mask foreground extent per straightened column, and the traced-path row of the straightened band.
var maskFirstYPx, maskLastYPx, maskMidYPx;
// 1 where the core mask genuinely reaches the column, 0 where the extent above was gap-filled.
var maskCorePresent;
// Output directory of the run currently in progress, so cleanup on an aborted run
// (from any helper function) can close windows and remove the empty folder.
var activeOutDir = "";

// Preference helpers to persist settings between runs
function prefGetNumber(key, def){ s = call("ij.Prefs.get", key, ""+def); return parseFloat(s); }
function prefGetBoolean(key, def){ if (def) defStr="true"; else defStr="false"; s = call("ij.Prefs.get", key, defStr); if (s == "true" || s == "1") return 1; return 0; }
function prefGetString(key, def){ return call("ij.Prefs.get", key, def); }
function prefSet(key, val){ call("ij.Prefs.set", key, ""+val); }

// Convert calibration units to um per pixel
function umPerPixelFromCalibration(){
    getPixelSize(u, cw, ch, cd);
    ump = cw;
    low = toLowerCase(u);
    if (startsWith(low,"nm")) ump = cw/1000.0;
    else if (startsWith(low,"mm")) ump = cw*1000.0;
    else if (startsWith(low,"cm")) ump = cw*10000.0;
    else if (indexOf(low,"pixel")>=0) ump = 1.0;
    return ump;
}

// ImageJ image IDs are often negative; 0 means not found.
function findImageIdByTitle(title) {
    if (title == "") return 0;
    esc = replace(title, "\\", "\\\\");
    esc = replace(esc, "\"", "\\\"");
    js = "var found=0; var ids=Packages.ij.WindowManager.getIDList();";
    js = js + "if(ids!=null){for(var i=0;i<ids.length;i++){var imp=Packages.ij.WindowManager.getImage(ids[i]);";
    js = js + "if(imp!=null&&imp.getTitle().equals(\"" + esc + "\")){found=ids[i];break;}}}";
    js = js + "''+found;";
    idStr = eval("script", js);
    if (idStr == "0" || idStr == "" || idStr == "null" || idStr == "NaN") return 0;
    return parseFloat(idStr);
}

function keeperIdForTitle(title) {
    if (title == "") return 0;
    if (analysisCorrWin != "" && title == analysisCorrWin) return idAnalysisCorr;
    if (analysis2CorrWin != "" && title == analysis2CorrWin) return idAnalysis2Corr;
    if (title == "Mask_core") return idMaskCore;
    if (maskTitle == "Mask_core" && title == maskTitle) return idMaskCore;
    return 0;
}

function isKeeperTitle(title) {
    if (title == "") return 0;
    if (analysisCorrWin != "" && title == analysisCorrWin) return 1;
    if (analysis2CorrWin != "" && title == analysis2CorrWin) return 1;
    if (title == "Mask_core") return 1;
    if (maskTitle == "Mask_core" && title == maskTitle) return 1;
    return 0;
}

function idStillOpen(id) {
    if (id == 0) return 0;
    if (isOpen(id)) return 1;
    return 0;
}

// Give the current batch-hidden image a real window; stays in batch mode.
function pinCurrentImage() {
    setBatchMode("show");
    return getImageID();
}

function keeperHasWindow(id) {
    if (id == 0) return 0;
    w = eval("script", "var imp=Packages.ij.WindowManager.getImage(" + id + "); if(imp!=null&&imp.getWindow()!=null) 1; else 0;");
    if (w == "1" || w == "1.0") return 1;
    return 0;
}

function showKeeperImage(title, id) {
    if (keeperHasWindow(id)) return 1;
    if (idStillOpen(id)) {
        selectImage(id);
        setBatchMode("show");
        return 1;
    }
    if (title != "" && activateImage(title)) {
        setBatchMode("show");
        return 1;
    }
    return 0;
}

// True for displayed windows and batch-hidden images (isOpen misses the latter).
function imageIsOpen(title) {
    if (title == "") return 0;
    if (nImages > 0) {
        if (getTitle() == title) return 1;
    }
    kid = keeperIdForTitle(title);
    if (idStillOpen(kid)) return 1;
    list = getList("image.titles");
    if (list.length > 0) {
        for (i = 0; i < list.length; i++) {
            if (list[i] == title) return 1;
        }
    }
    if (findImageIdByTitle(title) != 0) return 1;
    return 0;
}

// Select by image ID so batch-hidden images work. Skip toFront during hidden compute.
function activateImage(title) {
    if (title == "") return 0;
    if (nImages > 0) {
        if (getTitle() == title) return 1;
    }
    id = keeperIdForTitle(title);
    if (id == 0 || !idStillOpen(id)) id = findImageIdByTitle(title);
    if (id == 0) return 0;
    selectImage(id);
    if (hiddenComputeOn == 0) {
        eval("script", "var imp=Packages.ij.WindowManager.getImage(" + id + "); if(imp!=null&&imp.getWindow()!=null) imp.getWindow().toFront();");
    }
    return 1;
}

// Duplicate the current image and rename it. The duplicate is already current; do not selectWindow.
function duplicateCurrentAs(newTitle) {
    srcTitle = getTitle();
    srcId = getImageID();
    run("Duplicate...", "title=" + newTitle);
    if (getImageID() == srcId) {
        cancelRunAndExit(activeOutDir, "Could not duplicate \"" + srcTitle + "\" to \"" + newTitle + "\".");
    }
    if (getTitle() != newTitle) rename(newTitle);
}

// Select an image by title or exit with a clear message (works for batch-hidden images).
function requireWindow(title) {
    if (!activateImage(title)) {
        cancelRunAndExit(activeOutDir, "No window with the title \"" + title + "\" found.\n" +
             "Expected image: " + title);
    }
}

function closeWindowIfOpen(title) {
    if (isKeeperTitle(title)) return;
    id = findImageIdByTitle(title);
    if (id == 0) return;
    selectImage(id);
    setOption("Changes", false);
    close();
}

function closeDetectionTempsByPrefix() {
    js = "var ids=Packages.ij.WindowManager.getIDList(); var out='';";
    js = js + "if(ids!=null){for(var i=0;i<ids.length;i++){var imp=Packages.ij.WindowManager.getImage(ids[i]);";
    js = js + "if(imp!=null){var t=String(imp.getTitle());";
    js = js + "if(t.indexOf('Regions_')===0||t.indexOf('Result of')===0||t.indexOf('Segmented Particles')===0||t.indexOf('Drawing of')===0||t.indexOf('Mask of')===0)";
    js = js + "{if(out.length>0)out+=','; out+=ids[i];}}}} out;";
    idStr = eval("script", js);
    if (idStr == "" || idStr == "null") return;
    parts = split(idStr, ",");
    for (pi = 0; pi < parts.length; pi++) {
        id = parseFloat(parts[pi]);
        if (id == 0) continue;
        if (!idStillOpen(id)) continue;
        selectImage(id);
        t = getTitle();
        if (isKeeperTitle(t)) continue;
        setOption("Changes", false);
        close();
    }
}

// Hide newly created compute windows.
function beginHiddenCompute() {
    hiddenComputeOn = 1;
    setBatchMode(true);
}

// Close detection temps, then show remaining images.
function endHiddenCompute() {
    closeRegionsTempWindows();
    closeWindowIfOpen("Mask_seed");
    closeDetectionTempsByPrefix();
    if (!keeperHasWindow(idAnalysisCorr)) showKeeperImage(analysisCorrWin, idAnalysisCorr);
    if (!keeperHasWindow(idMaskCore)) showKeeperImage("Mask_core", idMaskCore);
    if (analysis2CorrWin != "" && !keeperHasWindow(idAnalysis2Corr)) showKeeperImage(analysis2CorrWin, idAnalysis2Corr);
    hiddenComputeOn = 0;
    setBatchMode(false);
}

// One select per image: mask flags, then source pixels. Fills profile/neuropil/valid per column
// and scanSumIn/scanCntIn/scanSumOut/scanCntOut for whole-image sums.
function fillMaskedColumnProfiles(srcWin, maskWin, sw, sh, profile, neuropil, valid) {
    nPix = sw * sh;
    if (nPix < 1) nPix = 1;
    maskFl = newArray(nPix);
    requireWindow(maskWin);
    for (x = 0; x < sw; x++) {
        colOff = x * sh;
        for (y = 0; y < sh; y++) {
            if (getPixel(x, y) > 0) maskFl[colOff + y] = 1;
            else maskFl[colOff + y] = 0;
        }
    }
    requireWindow(srcWin);
    scanSumIn = 0;
    scanCntIn = 0;
    scanSumOut = 0;
    scanCntOut = 0;
    for (x = 0; x < sw; x++) {
        colOff = x * sh;
        sumAIS = 0;
        cntAIS = 0;
        sumNP = 0;
        cntNP = 0;
        for (y = 0; y < sh; y++) {
            pv = getPixel(x, y);
            if (maskFl[colOff + y] == 1) {
                sumAIS += pv;
                cntAIS++;
            } else {
                sumNP += pv;
                cntNP++;
            }
        }
        if (cntAIS > 0) {
            profile[x] = sumAIS / cntAIS;
            valid[x] = 1;
        } else {
            profile[x] = 0;
            valid[x] = 0;
        }
        if (cntNP > 0) neuropil[x] = sumNP / cntNP;
        else neuropil[x] = 0;
        scanSumIn += sumAIS;
        scanCntIn += cntAIS;
        scanSumOut += sumNP;
        scanCntOut += cntNP;
    }
}

// Close all image windows whose title starts with prefix (e.g. peak preview after zoom adds "(800%)").
function closeWindowsWithTitlePrefix(prefix) {
    if (prefix == "") return;
    esc = replace(prefix, "\\", "\\\\");
    esc = replace(esc, "\"", "\\\"");
    js = "var ids=Packages.ij.WindowManager.getIDList(); var out='';";
    js = js + "if(ids!=null){for(var i=0;i<ids.length;i++){var imp=Packages.ij.WindowManager.getImage(ids[i]);";
    js = js + "if(imp!=null){var t=String(imp.getTitle());";
    js = js + "if(t.indexOf(\"" + esc + "\")===0){if(out.length>0)out+=','; out+=ids[i];}}}} out;";
    idStr = eval("script", js);
    if (idStr == "" || idStr == "null") return;
    parts = split(idStr, ",");
    for (pi = 0; pi < parts.length; pi++) {
        id = parseFloat(parts[pi]);
        if (id == 0) continue;
        if (!idStillOpen(id)) continue;
        selectImage(id);
        t = getTitle();
        if (isKeeperTitle(t)) continue;
        setOption("Changes", false);
        close();
    }
}

// saveAs often renames the window to the saved filename; close by expected title and basename.
function closeWindowAfterSave(expectedTitle, savePath) {
    if (expectedTitle != "") closeWindowIfOpen(expectedTitle);
    savedName = File.getName(savePath);
    dot = indexOf(savedName, ".");
    if (dot != -1) savedBase = substring(savedName, 0, dot); else savedBase = savedName;
    closeWindowsWithTitlePrefix(savedBase);
}

// ----- Cancellation cleanup -----
// ImageJ's macro language has no try/catch/finally, and clicking the built-in
// "Cancel" button on a Dialog (or waitForUser) aborts the macro immediately with
// no way to run code afterwards. To let the user abandon a poor-quality image and
// still have every window closed, we provide an explicit, catchable "Cancel run"
// action (see the tuning dialogs) plus these helpers, which the graceful exits use.

// Close every open image window without triggering "save changes?" prompts.
function closeAllImagesNoPrompt() {
    while (nImages > 0) {
        selectImage(nImages);
        setOption("Changes", false);
        close();
    }
}

// Delete the run's output folder (and its empty Analysis/Parameters/ROI subfolders)
// only if nothing has been written yet, so an abandoned run leaves no empty folders
// behind (which would otherwise clutter the auto-increment sequence).
function removeOutputDirIfEmpty(outDirToRemove) {
    if (outDirToRemove == "" || !File.exists(outDirToRemove)) return;
    subs = newArray("Analysis", "Parameters", "ROI");
    for (i = 0; i < subs.length; i++) {
        sp = outDirToRemove + File.separator + subs[i];
        if (File.exists(sp)) {
            l = getFileList(sp);
            if (l.length == 0) File.delete(sp);
        }
    }
    l2 = getFileList(outDirToRemove);
    if (l2.length == 0) File.delete(outDirToRemove);
}

// Close all analysis windows and remove the (empty) output folder for an aborted run.
function abortRunCleanup(outDirToRemove) {
    if (isOpen("ROI Manager")) { roiManager("reset"); selectWindow("ROI Manager"); run("Close"); }
    if (isOpen("Results")) { selectWindow("Results"); run("Close"); }
    closeAllImagesNoPrompt();
    hiddenComputeOn = 0;
    setBatchMode(false);
    removeOutputDirIfEmpty(outDirToRemove);
}

// Clean up and stop the macro. Pass msg="" for a silent, user-initiated cancel.
function cancelRunAndExit(outDirToRemove, msg) {
    abortRunCleanup(outDirToRemove);
    if (msg != "") exit(msg);
    showStatus("Run cancelled; output windows closed.");
    exit;
}

function backslashToSlash(p) {
    out = "";
    for (i = 0; i < lengthOf(p); i++) {
        c = substring(p, i, i + 1);
        if (c == "\\") out = out + "/";
        else out = out + c;
    }
    return out;
}

// Create dirPath and any missing parents (recursive). Exits with a clear message on failure.
function ensureDirectoryExists(dirPath) {
    if (dirPath == "") return;
    if (File.exists(dirPath)) return;
    javaPath = backslashToSlash(dirPath);
    eval("script", "new java.io.File(\"" + javaPath + "\").mkdirs();");
    if (!File.exists(dirPath)) {
        exit("Could not create output directory:\n" + dirPath +
             "\nCheck that the parent folder exists and is writable.");
    }
}

function ensureOutputSubdir(parentDir, subName) {
    ensureDirectoryExists(parentDir);
    path = parentDir + File.separator + subName;
    ensureDirectoryExists(path);
    return path;
}

function selectWindowByTitlePrefix(prefix) {
    if (isOpen(prefix)) {
        selectWindow(prefix);
        return;
    }
    titles = getList("image.titles");
    for (i = 0; i < titles.length; i++) {
        t = titles[i];
        if (startsWith(t, prefix)) {
            selectWindow(t);
            return;
        }
    }
}

// Next numbered output folder: imageStem_output_001, _002, ... (scans disk; resets if all deleted).
function allocateNumberedOutputDir(parentDir, imageStem) {
    prefix = imageStem + "_output_";
    maxRun = 0;
    list = getFileList(parentDir);
    for (i = 0; i < list.length; i++) {
        name = list[i];
        if (!endsWith(name, "/") && !endsWith(name, "\\")) continue;
        folderName = name;
        if (endsWith(folderName, "/")) folderName = substring(folderName, 0, lengthOf(folderName) - 1);
        else if (endsWith(folderName, "\\")) folderName = substring(folderName, 0, lengthOf(folderName) - 1);
        if (startsWith(folderName, prefix)) {
            suffix = substring(folderName, lengthOf(prefix), lengthOf(folderName));
            if (lengthOf(suffix) == 3) {
                runN = parseInt(suffix);
                if (runN > maxRun) maxRun = runN;
            }
        }
    }
    nextRun = maxRun + 1;
    runTag = "" + nextRun;
    while (lengthOf(runTag) < 3) runTag = "0" + runTag;
    outPath = parentDir + imageStem + "_output_" + runTag;
    ensureDirectoryExists(outPath);
    return outPath;
}

// Replace characters that are illegal in Windows/macOS/Linux file or folder names.
function sanitizeForFilename(s) {
    out = "";
    for (i = 0; i < lengthOf(s); i++) {
        c = substring(s, i, i + 1);
        if (c == "\\" || c == "/" || c == ":" || c == "*" || c == "?" || c == "\"" || c == "<" || c == ">" || c == "|") out = out + "_";
        else out = out + c;
    }
    while (endsWith(out, " ") || endsWith(out, ".")) out = substring(out, 0, lengthOf(out) - 1);
    while (startsWith(out, " ")) out = substring(out, 1, lengthOf(out));
    return out;
}

// Output base name for a run.
//  * Ordinary single-image files: the file name without its extension (unchanged
//    behaviour, using the first dot exactly as before).
//  * Multi-image container files (e.g. Leica .lif): Bio-Formats titles the opened
//    window "<file>.lif - <series name>", so we adopt the individual series/stack
//    name instead of the shared container name. Single-image files never carry this
//    " - <series>" suffix, so their naming is completely unaffected.
function deriveOutputBaseName(inPath, windowTitle) {
    fileName = File.getName(inPath);
    fileBase = fileName;
    dot = indexOf(fileName, ".");
    if (dot != -1) fileBase = substring(fileName, 0, dot);
    sepFull = fileName + " - ";
    sepBase = fileBase + " - ";
    series = "";
    if (startsWith(windowTitle, sepFull)) series = substring(windowTitle, lengthOf(sepFull), lengthOf(windowTitle));
    else if (startsWith(windowTitle, sepBase)) series = substring(windowTitle, lengthOf(sepBase), lengthOf(windowTitle));
    series = sanitizeForFilename(series);
    if (series != "") return series;
    return fileBase;
}

function clearImageRoi() {
    run("Select None");
    eval("script", "var imp=Packages.ij.IJ.getImage(); if(imp!=null) imp.deleteRoi();");
}

// Checkbox + right-aligned hint for Output Metrics dialog rows.
function addMetricCheckbox(label, defaultVal, hint) {
    padW = 42;
    padded = label;
    while (lengthOf(padded) < padW) padded = padded + " ";
    Dialog.addCheckbox(padded, defaultVal);
    Dialog.addToSameRow();
    Dialog.addMessage(hint);
}

// Save plain straightened PNGs before closing intermediate windows (uses macro outDir, base).
function saveStraightenedPlainPngs() {
    if (imageIsOpen(roiStraightWin)) {
        activateImage(roiStraightWin);
        run("Remove Overlay");
        duplicateCurrentAs(roiPlainWin);
        run("Flatten");
        ch1PlainPath = roiOutDir + File.separator + base + roiExportTag + "_straight.png";
        saveAs("PNG", ch1PlainPath);
        closeWindowAfterSave(roiPlainWin, ch1PlainPath);
    }
    if (imageIsOpen(analysisStraightWin)) {
        activateImage(analysisStraightWin);
        run("Remove Overlay");
        duplicateCurrentAs(analysisPlainWin);
        run("Flatten");
        ch2PlainPath = analysisOutDir + File.separator + base + analysisExportTag + "_straight.png";
        saveAs("PNG", ch2PlainPath);
        closeWindowAfterSave(analysisPlainWin, ch2PlainPath);
    }
    if (adjDualMode && imageIsOpen(analysis2StraightWin)) {
        activateImage(analysis2StraightWin);
        run("Remove Overlay");
        duplicateCurrentAs(analysis2PlainWin);
        run("Flatten");
        ch2bPlainPath = analysisOutDir + File.separator + base + analysis2ExportTag + "_straight.png";
        saveAs("PNG", ch2bPlainPath);
        closeWindowAfterSave(analysis2PlainWin, ch2bPlainPath);
    }
}

// Close windows no longer needed once profiles/CV are in arrays (keep Mask_core / maskTitle, analysisCorrWin).
function closePrePeakTuningWindows() {
    if (imageIsOpen(analysisCorrWin)) {
        activateImage(analysisCorrWin);
        run("Select None");
    }
    if (adjDualMode && imageIsOpen(analysis2CorrWin)) {
        activateImage(analysis2CorrWin);
        run("Select None");
    }
    closeWindowIfOpen(roiStraightWin);
    closeWindowIfOpen(roiCorrWin);
    closeWindowIfOpen(analysisStraightWin);
    if (analysis2StraightWin != "") closeWindowIfOpen(analysis2StraightWin);
    closeWindowIfOpen(analysisNormWin);
    closeWindowIfOpen("Mask_straight");
    closeWindowIfOpen("ArcProfile_FM");
}

// Build binary Mask_work from Mask_straight using threshold fraction thrFrac (uses rbUm, maskDilUm, umPerPx).
function buildMaskWorkFromStraight(thrFrac) {
    if (isOpen("Mask_work")) { selectWindow("Mask_work"); run("Close"); }
    selectWindow("Mask_straight");
    run("Duplicate...", "title=Mask_work");
    selectWindow("Mask_work");
    if (rbUm > 0) {
        rbPxM = round(rbUm / umPerPx);
        if (rbPxM < 1) rbPxM = 1;
        run("Subtract Background...", "rolling=" + rbPxM + " sliding");
    }
    makeRectangle(0, 0, getWidth(), getHeight());
    getStatistics(aMask, meanMask, minMask, maxMask, stdMask);
    if (maxMask <= 0) maxMask = 1;
    thrVal = thrFrac * maxMask;
    setOption("BlackBackground", true);
    setThreshold(thrVal, maxMask);
    run("Convert to Mask");
    if (maskDilUm > 0) {
        dil_px = round(maskDilUm / umPerPx);
        if (dil_px > 0) run("Dilate", "iterations=" + dil_px);
    }
}

// Assign split channel to role: rename on first use of channel index, duplicate if index reused.
function assignRoleChannel(chIndex, srcTitle, roleTitle, consumedTitles) {
    if (consumedTitles[chIndex] == "") {
        requireWindow(srcTitle);
        rename(roleTitle);
        consumedTitles[chIndex] = roleTitle;
    } else {
        requireWindow(consumedTitles[chIndex]);
        // "duplicate" keeps every slice; without it ImageJ copies only the current plane (breaks z-stacks).
        if (nSlices > 1) run("Duplicate...", "title=" + roleTitle + " duplicate");
        else run("Duplicate...", "title=" + roleTitle);
    }
}

// Close unsplit stack and leftover C{n}- split windows after role assignment.
function closeLeftoverChannelWindows(origTitle, nChannels) {
    if (isOpen(origTitle)) {
        selectWindow(origTitle);
        run("Close");
    }
    titlesAll = getList("image.titles");
    for (ci = 1; ci <= nChannels; ci++) {
        cPrefix = "C" + d2s(ci, 0) + "-";
        for (i = 0; i < titlesAll.length; i++) {
            t = titlesAll[i];
            if (startsWith(t, cPrefix) && isOpen(t)) {
                selectWindow(t);
                run("Close");
            }
        }
    }
}

// Exact split-channel title for 1-based index (C{n}-..., avoids C1 matching C10).
function findSplitChannelTitle(chIdx) {
    prefix = "C" + d2s(chIdx, 0) + "-";
    titlesList = getList("image.titles");
    for (i = 0; i < titlesList.length; i++) {
        t = titlesList[i];
        if (startsWith(t, prefix)) return t;
    }
    return "";
}

// Draw mask outline on a duplicate of ch2 for non-blocking preview (no interactive Threshold).
// Uses ROI Manager briefly so the selection transfers between mask and preview windows.
function buildMaskPreview(ch2Title, maskTitle, previewTitle) {
    if (isOpen(previewTitle)) { selectWindow(previewTitle); run("Close"); }
    selectWindow(ch2Title);
    run("Duplicate...", "title=" + previewTitle);
    selectWindow(previewTitle);
    run("RGB Color");
    selectWindow(maskTitle);
    setOption("BlackBackground", true);
    setThreshold(1, 255);
    run("Create Selection");
    if (selectionType() == -1) {
        showMessage("Mask preview", "Mask is empty at this threshold.\nTry a lower threshold fraction, then choose Update preview.");
        selectWindow(previewTitle);
        run("Select None");
        return;
    }
    roiManager("Add");
    maskRoiIdx = roiManager("count") - 1;
    selectWindow(previewTitle);
    roiManager("Select", maskRoiIdx);
    setLineWidth(2);
    setForegroundColor(255, 0, 0);
    run("Draw");
    roiManager("Select", maskRoiIdx);
    roiManager("Delete");
    run("Select None");
}

function tuckStraightSourceWindows() {
    tuckX = -5000;
    tuckY = -5000;
    if (isOpen("Mask_straight")) {
        selectWindow("Mask_straight");
        setLocation(tuckX, tuckY);
    }
    if (isOpen(analysisStraightWin)) {
        selectWindow(analysisStraightWin);
        setLocation(tuckX, tuckY);
    }
    if (adjDualMode && analysis2StraightWin != "" && isOpen(analysis2StraightWin)) {
        selectWindow(analysis2StraightWin);
        setLocation(tuckX, tuckY);
    }
}

// Side-by-side tuning previews; mask-channel window on top. Dialog is still modal.
function layoutTuningPreviews() {
    if (!isOpen(maskPreviewWin) || !isOpen(analysisMaskPreviewWin)) return;
    padX = 40;
    padY = 80;
    selectWindow(maskPreviewWin);
    setLocation(padX, padY);
    mw = getWidth();
    mh = getHeight();
    selectWindow(analysisMaskPreviewWin);
    setLocation(padX + mw + 20, padY);
    if (adjDualMode && analysis2MaskPreviewWin != "" && isOpen(analysis2MaskPreviewWin)) {
        selectWindow(analysis2MaskPreviewWin);
        setLocation(padX + 2 * (mw + 20), padY);
    }
    // selectWindow activates and raises the window (no run("Bring to Front") in macros).
    tuckStraightSourceWindows();
    selectWindow(maskPreviewWin);
}

function syncBgsubPreviewDisplay(rawWin, previewWin) {
    selectWindow(rawWin);
    getStatistics(area, mean, minRaw, maxRaw, std);
    setMinAndMax(minRaw, maxRaw);
    if (isOpen(previewWin)) {
        selectWindow(previewWin);
        setMinAndMax(minRaw, maxRaw);
    }
}

function buildBgsubPreview(srcWin, previewTitle, rbUmVal) {
    closeWindowIfOpen(previewTitle);
    selectWindow(srcWin);
    run("Duplicate...", "title=" + previewTitle);
    selectWindow(previewTitle);
    if (rbUmVal > 0) {
        rbPx = round(rbUmVal / umPerPx);
        if (rbPx < 1) rbPx = 1;
        run("Subtract Background...", "rolling=" + rbPx + " sliding");
    }
    syncBgsubPreviewDisplay(srcWin, previewTitle);
}

function layoutBgsubPreviews() {
    if (!isOpen(roiBgsubPreviewWin) || !isOpen(analysisBgsubPreviewWin)) return;
    padX = 40;
    padY = 80;
    selectWindow(roiBgsubPreviewWin);
    setLocation(padX, padY);
    mw = getWidth();
    selectWindow(analysisBgsubPreviewWin);
    setLocation(padX + mw + 20, padY);
    syncBgsubPreviewDisplay("Mask_straight", roiBgsubPreviewWin);
    syncBgsubPreviewDisplay(analysisStraightWin, analysisBgsubPreviewWin);
    if (adjDualMode && analysis2BgsubPreviewWin != "" && isOpen(analysis2BgsubPreviewWin)) {
        selectWindow(analysis2BgsubPreviewWin);
        setLocation(padX + 2 * (mw + 20), padY);
        syncBgsubPreviewDisplay(analysis2StraightWin, analysis2BgsubPreviewWin);
    }
    tuckStraightSourceWindows();
    selectWindow(roiBgsubPreviewWin);
}

// Per-column max of ch2 inside mask foreground. Select mask once, then analysis once.
function buildArcProfileFromMask(maskTitle, ch2Win, swArc, shArc, profileArc, yAtMax, arcValid) {
    nPix = swArc * shArc;
    if (nPix < 1) nPix = 1;
    maskFl = newArray(nPix);
    requireWindow(maskTitle);
    for (x = 0; x < swArc; x++) {
        colOff = x * shArc;
        for (y = 0; y < shArc; y++) {
            if (getPixel(x, y) > 0) maskFl[colOff + y] = 1;
            else maskFl[colOff + y] = 0;
        }
    }
    requireWindow(ch2Win);
    for (x = 0; x < swArc; x++) {
        colOff = x * shArc;
        maxV = -1e10;
        bestY = 0;
        found = 0;
        for (y = 0; y < shArc; y++) {
            if (maskFl[colOff + y] == 1) {
                pv = getPixel(x, y);
                if (pv > maxV) {
                    maxV = pv;
                    bestY = y;
                    found = 1;
                }
            }
        }
        if (found) {
            profileArc[x] = maxV;
            yAtMax[x] = bestY;
            arcValid[x] = 1;
        } else {
            profileArc[x] = 0;
            yAtMax[x] = 0;
            arcValid[x] = 0;
        }
    }
}

// Refine peak positions using mask + Ch2 pixels in +/- refinePx columns.
function refinePeaksInMask(xs, ys, vals, nPeaks, maskTitle, ch2Win, swArc, shArc, refinePx) {
    if (nPeaks < 1) return;
    maskCol = newArray(shArc);
    for (i = 0; i < nPeaks; i++) {
        bestV = vals[i];
        bestX = xs[i];
        bestY = ys[i];
        x0 = xs[i] - refinePx;
        if (x0 < 0) x0 = 0;
        x1 = xs[i] + refinePx;
        if (x1 >= swArc) x1 = swArc - 1;
        for (xr = x0; xr <= x1; xr++) {
            requireWindow(maskTitle);
            for (y = 0; y < shArc; y++) {
                if (getPixel(xr, y) > 0) maskCol[y] = 1; else maskCol[y] = 0;
            }
            requireWindow(ch2Win);
            for (y = 0; y < shArc; y++) {
                if (maskCol[y]) {
                    pv = getPixel(xr, y);
                    if (pv > bestV) {
                        bestV = pv;
                        bestX = xr;
                        bestY = y;
                    }
                }
            }
        }
        xs[i] = bestX;
        ys[i] = bestY;
        vals[i] = bestV;
    }
}

function computeArcTolAbs(profileArc, arcValid, n, tolPct) {
    minArc = 0;
    maxArc = 0;
    arcHasValid = 0;
    for (x = 0; x < n; x++) {
        if (arcValid[x]) {
            maxV = profileArc[x];
            if (arcHasValid == 0) {
                minArc = maxV;
                maxArc = maxV;
                arcHasValid = 1;
            } else {
                if (maxV < minArc) minArc = maxV;
                if (maxV > maxArc) maxArc = maxV;
            }
        }
    }
    if (arcHasValid == 0) {
        minArc = 0;
        maxArc = 1;
    }
    rangeArc = maxArc - minArc;
    if (rangeArc <= 0) rangeArc = 1;
    return (tolPct / 100.0) * rangeArc;
}

function findPeaksArcFiji(profileArc, arcValid, n, tolAbs, peaksBright) {
    peakX = newArray(0);
    if (n < 1) return peakX;
    closeWindowIfOpen("ArcProfile_FM");
    newImage("ArcProfile_FM", "32-bit", n, 1, 1);
    selectWindow("ArcProfile_FM");
    for (x = 0; x < n; x++) {
        if (arcValid[x]) setPixel(x, 0, profileArc[x]);
        else setPixel(x, 0, 0);
    }
    if (!peaksBright) run("Invert");
    if (isOpen("Results")) { selectWindow("Results"); run("Clear Results"); }
    run("Find Maxima...", "noise=" + d2s(tolAbs, 6) + " light output=[List]");
    nFmRaw = nResults;
    nOut = 0;
    peakX = newArray(nFmRaw);
    for (i = 0; i < nFmRaw; i++) {
        px = round(getResult("X", i));
        if (px < 0) px = 0;
        if (px >= n) px = n - 1;
        if (arcValid[px]) {
            peakX[nOut] = px;
            nOut++;
        }
    }
    closeWindowIfOpen("ArcProfile_FM");
    if (isOpen("Results")) { selectWindow("Results"); run("Close"); }
    if (nOut == 0) return newArray(0);
    if (nOut == nFmRaw) return peakX;
    out = newArray(nOut);
    for (k = 0; k < nOut; k++) out[k] = peakX[k];
    return out;
}

// Fill rangeOut[0]=min, rangeOut[1]=max, rangeOut[2]=range for Ch2 inside mask (selection-based, fast).
function computeMaskedCh2Range(maskTitle, ch2Win, rangeOut) {
    closeWindowIfOpen("Range_proc");
    requireWindow(ch2Win);
    duplicateCurrentAs("Range_proc");
    imageCalculator("Multiply create 32-bit", "Range_proc", maskTitle);
    if (imageIsOpen("Result of Range_proc")) {
        activateImage("Result of Range_proc");
        rename("Range_proc");
    }
    requireWindow(maskTitle);
    setOption("BlackBackground", true);
    setThreshold(1, 255);
    run("Create Selection");
    if (selectionType() == -1) {
        closeWindowIfOpen("Range_proc");
        rangeOut[0] = 0;
        rangeOut[1] = 1;
        rangeOut[2] = 1;
        return;
    }
    roiManager("reset");
    roiManager("Add");
    activateImage("Range_proc");
    roiManager("Select", 0);
    getStatistics(area, mean, minV, maxV, std);
    roiManager("reset");
    run("Select None");
    closeWindowIfOpen("Range_proc");
    rangeV = maxV - minV;
    if (rangeV <= 0) rangeV = 1;
    rangeOut[0] = minV;
    rangeOut[1] = maxV;
    rangeOut[2] = rangeV;
}

function isArcPeakMethod(method) {
    if (indexOf(method, "Arc profile") >= 0) return 1;
    return 0;
}

function sumRegionAreas(nReg) {
    sumA = 0;
    for (i = 0; i < nReg; i++) sumA += peakRegAreaWork[i];
    return sumA;
}

// Absolute intensity threshold from in-mask range and tolerance percent.
function thresholdAbsFromPct(minV, maxV, tolPct, peaksBright) {
    rangeV = maxV - minV;
    if (rangeV <= 0) rangeV = 1;
    if (peaksBright) return minV + (tolPct / 100.0) * rangeV;
    return maxV - (tolPct / 100.0) * rangeV;
}

// Baseline-relative threshold for the flattened residual: tol% of the way from the
// in-mask mean (baseline ~ 0 after flattening) toward the brightest (or darkest) patch.
function thresholdAbsFromMeanMax(meanV, minV, maxV, tolPct, peaksBright) {
    if (peaksBright) {
        span = maxV - meanV;
        if (span <= 0) span = 1;
        return meanV + (tolPct / 100.0) * span;
    }
    span = meanV - minV;
    if (span <= 0) span = 1;
    return meanV - (tolPct / 100.0) * span;
}

// Zero negative pixels after baseline flatten; keep only fluorescence above local baseline.
function clampPositiveResidual(winTitle) {
    requireWindow(winTitle);
    getDimensions(cw, ch, cc, cs, cf);
    for (y = 0; y < ch; y++) {
        for (x = 0; x < cw; x++) {
            v = getPixel(x, y);
            if (v < 0) setPixel(x, y, 0);
        }
    }
}

// Copy spatial calibration from one image window to another (for consistent BX/BY units).
// Uses getVoxelSize (valid macro function); getPixelWidth/getUnit are NOT macro functions.
function syncImageCalibration(sourceTitle, targetTitle) {
    requireWindow(sourceTitle);
    getVoxelSize(calPw, calPh, calVd, calUnit);
    requireWindow(targetTitle);
    run("Properties...", "channels=1 slices=1 frames=1 pixel_width=" + calPw + " pixel_height=" + calPh + " voxel_depth=" + calVd + " unit=" + calUnit);
}

// 32-bit duplicate of analysis channel for particle intensity redirect (keeps calibration).
function prepareRegionsMeasureImage(ch2Win, measureTitle) {
    closeWindowIfOpen(measureTitle);
    requireWindow(ch2Win);
    run("Select None");
    duplicateCurrentAs(measureTitle);
    if (bitDepth() != 32) run("32-bit");
}

// Force a window's spatial calibration to umPerPx so Analyze Particles reports geometry
// (BX/BY/Width/Height/Area/Feret) in microns regardless of whether the source used image
// calibration or a manual FOV width. Only affects reported units, not detection.
function setMicronScale(winTitle, umPerPxVal) {
    if (!imageIsOpen(winTitle)) return;
    activateImage(winTitle);
    // setVoxelSize only ever affects the current image; unlike Set Scale... it has no global option.
    setVoxelSize(umPerPxVal, umPerPxVal, 1, "um");
}

// Restore the region Set Measurements config with no redirect (call after region particle analysis).
function resetRegionMeasurementsRedirect() {
    run("Set Measurements...", "area mean bounding feret's shape redirect=None decimal=3");
}

// Remap label/watershed mask (values 1..N) to true binary 0/255 for Analyze Particles.
function binarizeNonzeroPixels(winTitle) {
    requireWindow(winTitle);
    getDimensions(bw, bh, bc, bs, bf);
    for (y = 0; y < bh; y++) {
        for (x = 0; x < bw; x++) {
            if (getPixel(x, y) > 0) setPixel(x, y, 255);
        }
    }
}

// Count particles on a binary mask without storing centroids (for tuning status).
function countParticlesOnBinary(maskTitle, ch2Win, minAreaUm2, maxAreaUm2, umPerPx) {
    requireWindow(maskTitle);
    um2PerPx = umPerPx * umPerPx;
    minPx = round(minAreaUm2 / um2PerPx);
    if (minPx < 1) minPx = 1;
    if (maxAreaUm2 > 0) maxPx = round(maxAreaUm2 / um2PerPx);
    else maxPx = 99999999;
    sizeStr = "" + minPx + "-" + maxPx;
    if (isOpen("Results")) { selectWindow("Results"); run("Clear Results"); }
    activateImage(maskTitle);
    run("Analyze Particles...", "size=" + sizeStr + " pixel circularity=0.00-1.00 show=Nothing display clear redirect=" + ch2Win);
    nCnt = nResults;
    if (isOpen("Results")) { selectWindow("Results"); run("Close"); }
    if (imageIsOpen(maskTitle)) { activateImage(maskTitle); run("Select None"); }
    if (imageIsOpen(ch2Win)) { activateImage(ch2Win); run("Select None"); }
    return nCnt;
}

// Count all connected components on a binary mask (no area filter; for tuning diagnostics).
function countAllParticlesOnBinary(maskTitle, ch2Win) {
    requireWindow(maskTitle);
    if (isOpen("Results")) { selectWindow("Results"); run("Clear Results"); }
    activateImage(maskTitle);
    run("Analyze Particles...", "size=0-Infinity circularity=0.00-1.00 show=Nothing display clear redirect=" + ch2Win);
    nCnt = nResults;
    if (isOpen("Results")) { selectWindow("Results"); run("Close"); }
    if (imageIsOpen(maskTitle)) { activateImage(maskTitle); run("Select None"); }
    if (imageIsOpen(ch2Win)) { activateImage(ch2Win); run("Select None"); }
    return nCnt;
}

// Store bounding boxes of all components (no area filter) for diagnostic preview, and
// persist the true particle outlines to a temp ROI set so the preview can trace them.
function captureAllRegionsFromBinary(maskTitle, ch2Win, regBX, regBY, regW, regH, maxN) {
    requireWindow(maskTitle);
    if (isOpen("Results")) { selectWindow("Results"); run("Clear Results"); }
    roiManager("reset");
    regionDiagRoiPath = "";
    run("Set Measurements...", "area bounding redirect=None decimal=3");
    activateImage(maskTitle);
    run("Analyze Particles...", "size=0-Infinity circularity=0.00-1.00 show=Nothing display clear add redirect=" + ch2Win);
    nCap = nResults;
    if (nCap > maxN) nCap = maxN;
    for (i = 0; i < nCap; i++) {
        roiManager("Select", i);
        getSelectionBounds(regBX[i], regBY[i], regW[i], regH[i]);
    }
    // Save the speckle outlines (still in the ROI Manager) for the diagnostic preview.
    if (nCap > 0) {
        diagPath = getDirectory("temp") + "neurite_regionDiag.zip";
        roiManager("Deselect");
        roiManager("Save", diagPath);
        regionDiagRoiPath = diagPath;
    }
    if (imageIsOpen(maskTitle)) { activateImage(maskTitle); run("Select None"); }
    if (imageIsOpen(ch2Win)) { activateImage(ch2Win); run("Select None"); }
    if (isOpen("Results")) { selectWindow("Results"); run("Close"); }
    return nCap;
}

// Save pixel-accurate outlines of kept regions for final overlay export.
function saveKeptRegionExportRois(nRegions) {
    if (nRegions < 1) return;
    nAll = roiManager("count");
    if (nAll < 1) return;
    tempAllPath = getDirectory("temp") + "neurite_regionAll.zip";
    roiManager("Deselect");
    roiManager("Save", tempAllPath);
    roiManager("reset");
    roiManager("Open", tempAllPath);
    nAll = roiManager("count");
    for (j = nAll - 1; j >= 0; j--) {
        keep = 0;
        for (i = 0; i < nRegions; i++) {
            if (peakRegionKeptRoiIdx[i] == j) keep = 1;
        }
        if (keep == 0) {
            roiManager("Select", j);
            roiManager("Delete");
        }
    }
    if (roiManager("count") > 0) {
        exportPath = getDirectory("temp") + "neurite_regionExport.zip";
        roiManager("Save", exportPath);
        regionExportRoiPath = exportPath;
    }
    roiManager("reset");
}

// After detection is accepted (interactive) or auto-run (batch): save region ROI outlines and arc profile export arrays.
function finalizePeakDetectionResults() {
    if (isRegionMethodUsed && nPeaks > 0) {
        if (regionExportRoiPath == "" || !File.exists(regionExportRoiPath)) {
            saveKeptRegionExportRois(nPeaks);
        }
    }
    buildArcProfileFromMask(maskTitle, analysisCorrWin, swArc, shArc, profileArc, yAtMax, arcValid);
    profileArcExport = newArray(swArc);
    arcValidExport = newArray(swArc);
    swArcExport = swArc;
    arcProfMin = 1e30;
    arcProfMax = -1e30;
    arcProfValidCols = 0;
    for (x = 0; x < swArc; x++) {
        profileArcExport[x] = profileArc[x];
        arcValidExport[x] = arcValid[x];
        if (arcValid[x]) {
            arcProfValidCols++;
            if (profileArc[x] < arcProfMin) arcProfMin = profileArc[x];
            if (profileArc[x] > arcProfMax) arcProfMax = profileArc[x];
        }
    }
    if (arcProfValidCols == 0) {
        arcProfMin = 0;
        arcProfMax = 0;
    }
}

// Transfer the neurite mask boundary to targetTitle as its current selection. The core mask is
// used rather than the full one, so the outline shows exactly what the distances are measured from.
// Returns 1 when a selection was restored, 0 when the mask has no foreground.
function restoreMaskOutlineSelection(targetTitle) {
    outlineSrc = adjCoreMaskTitle;
    if (outlineSrc == "" || !imageIsOpen(outlineSrc)) outlineSrc = maskTitle;
    if (!imageIsOpen(outlineSrc) || !imageIsOpen(targetTitle)) return 0;
    activateImage(outlineSrc);
    setOption("BlackBackground", true);
    setThreshold(1, 255);
    run("Create Selection");
    hasSel = (selectionType() != -1);
    resetThreshold();
    run("Select None");
    if (!hasSel) return 0;
    activateImage(targetTitle);
    run("Restore Selection");
    return 1;
}

// Adjacent-particle mode: trace the mask boundary in magenta so the distance of each
// particle from the neurite can be judged by eye. No-op in every other mode.
function drawMaskOutlineOnWindow(targetTitle) {
    if (!adjacentMode) return;
    if (restoreMaskOutlineSelection(targetTitle) == 0) return;
    setForegroundColor(255, 0, 255);
    setLineWidth(1);
    run("Draw");
    run("Select None");
}

function addMaskOutlineOverlay(targetTitle) {
    if (!adjacentMode) return;
    if (restoreMaskOutlineSelection(targetTitle) == 0) return;
    run("Add Selection...", "stroke=magenta width=1");
    run("Select None");
}

// Draw a single horizontal guide line across the image, skipping rows outside it.
function drawGuideRow(imgW, imgH, row) {
    if (row < 0 || row > imgH - 1) return;
    makeLine(0, row, imgW - 1, row);
    run("Draw");
    run("Select None");
}

// Adjacent-particle mode: mark the straightened midline in green along with the two edges of the
// acceptance band, so the rejection threshold is visible against the particles while tuning.
// No-op in every other mode.
function drawAdjacentGuideLinesOnWindow(targetTitle) {
    if (!adjacentMode) return;
    if (!imageIsOpen(targetTitle)) return;
    activateImage(targetTitle);
    run("Select None");
    getDimensions(gW, gH, gC, gS, gF);
    setForegroundColor(0, 255, 0);
    setLineWidth(1);
    midRow = round(maskMidYPx);
    drawGuideRow(gW, gH, midRow);
    if (adjUmPerPx > 0) {
        bandPx = round(adjMaxDistMidUm / adjUmPerPx);
        if (bandPx > 0) {
            drawGuideRow(gW, gH, midRow - bandPx);
            drawGuideRow(gW, gH, midRow + bandPx);
        }
    }
    run("Select None");
}

// Add peak/region overlays to the active image window (used for final PNG export).
function addPeakOverlaysToImage(isRegions, nPeaks, xs, ys, markerUm, umPerPx) {
    overlayTarget = getTitle();
    addMaskOutlineOverlay(overlayTarget);
    activateImage(overlayTarget);
    if (nPeaks < 1) return;
    // Add Selection remembers the stroke from its last use, so state it explicitly: without it the
    // patches inherit the magenta of the mask outline added just above. Adjacent mode draws them as
    // thin as that outline so the two read as separate features rather than one blob.
    if (adjacentMode && adjDualMode && adjPassUseAdj2Prefs) patchStroke = "stroke=cyan width=1";
    else if (adjacentMode) patchStroke = "stroke=yellow width=1";
    else patchStroke = "stroke=yellow width=3";
    if (isRegions) {
        drew = 0;
        if (regionExportRoiPath == "" && nPeaks > 0 && roiManager("count") > 0) {
            saveKeptRegionExportRois(nPeaks);
        }
        roiManager("reset");
        if (regionExportRoiPath != "" && File.exists(regionExportRoiPath)) {
            roiManager("Open", regionExportRoiPath);
            nExp = roiManager("count");
            for (i = 0; i < nExp; i++) {
                roiManager("Select", i);
                run("Add Selection...", patchStroke);
                run("Select None");
            }
            roiManager("reset");
            drew = 1;
        }
        if (drew == 0) {
            showStatus("Region overlay: no saved ROI outlines — skipping misaligned bounding boxes.");
        }
    } else {
        radPx = round((markerUm / 2.0) / umPerPx);
        if (radPx < 1) radPx = 1;
        for (i = 0; i < nPeaks; i++) {
            cx = xs[i];
            cy = ys[i];
            makeOval(cx - radPx, cy - radPx, 2 * radPx, 2 * radPx);
            run("Add Selection...", patchStroke);
            run("Select None");
        }
    }
    if (Overlay.size > 0) Overlay.drawLabels(false);
}

function closeRegionsTempWindows() {
    closeWindowIfOpen("Regions_proc");
    closeWindowIfOpen("Regions_binary");
    closeWindowIfOpen("Regions_binary_pre");
    closeWindowIfOpen("Regions_intensity");
    closeWindowIfOpen("Regions_measure");
    closeWindowIfOpen("Regions_src");
    closeWindowIfOpen("Regions_small");
    closeWindowIfOpen("Regions_base");
    closeDetectionTempsByPrefix();
}

// Threshold + Analyze Particles; fills xs/ys/vals, peakReg*Work, and ROI Manager outlines.
// Adjacent mode never uses baseline flatten (absolute intensity on the full band).
function findRegions2DMasked(maskTitle, ch2Win, tolPct, peaksBright, minAreaUm2, maxAreaUm2, smoothSigmaUm, umPerPx, xs, ys, vals) {
    if (adjacentMode) regionDogUm = 0;
    closeRegionsTempWindows();
    roiManager("reset");
    regionCountThr = 0;
    regionCountThrFrag = 0;
    regionCountWs = 0;
    regionCountFrag = 0;
    regionResidualSpan = 0;
    regionCountPreDist = 0;
    regionDiagN = 0;

    // Build the working intensity image. Flatten the neurite's non-zero baseline so that
    // patches (locally brighter than the surrounding neurite) stand out as a positive residual.
    // Flattening runs on the full Ch2 image (not the masked one) to avoid dark-edge rim artifacts.
    isFlatten = (regionDogUm > 0);
    requireWindow(ch2Win);
    duplicateCurrentAs("Regions_src");
    if (isFlatten) {
        // High-pass residual = small-scale (lightly denoised) image minus large-scale baseline.
        sigmaSmall = 0;
        if (smoothSigmaUm > 0) sigmaSmall = smoothSigmaUm / umPerPx;
        // Baseline width (um): larger than patch size so patches survive in the residual.
        sigmaLarge = regionDogUm / umPerPx;
        if (sigmaLarge < 1.0) sigmaLarge = 1.0;
        activateImage("Regions_src");
        duplicateCurrentAs("Regions_small");
        if (sigmaSmall > 0) { run("32-bit"); run("Gaussian Blur...", "sigma=" + sigmaSmall); }
        else run("32-bit");
        activateImage("Regions_src");
        duplicateCurrentAs("Regions_base");
        run("32-bit");
        run("Gaussian Blur...", "sigma=" + sigmaLarge);
        imageCalculator("Subtract create 32-bit", "Regions_small", "Regions_base");
        closeWindowIfOpen("Regions_src");
        closeWindowIfOpen("Regions_small");
        closeWindowIfOpen("Regions_base");
        if (imageIsOpen("Result of Regions_small")) { activateImage("Result of Regions_small"); rename("Regions_proc"); }
        // Bright patches sit above the local baseline (positive residual); clamp the
        // negatives away so threshold/stats see only the above-baseline excess. For dark
        // patches the informative signal is the negative residual, so leave it unclamped.
        if (peaksBright) clampPositiveResidual("Regions_proc");
    } else {
        // Legacy raw-intensity path (no baseline flattening).
        activateImage("Regions_src");
        run("32-bit");
        if (smoothSigmaUm > 0) {
            sigmaPx = smoothSigmaUm / umPerPx;
            if (sigmaPx < 0.1) sigmaPx = 0.1;
            run("Gaussian Blur...", "sigma=" + sigmaPx);
        }
        rename("Regions_proc");
    }

    // Statistics of the working image that set the threshold. Normally restricted to the mask
    // (selection-based; no mask-multiply scaling). In adjacent-particle mode the signal of
    // interest lies outside the mask, so the whole straightened band is used instead.
    rangeOut = newArray(3);
    if (adjacentMode) {
        activateImage("Regions_proc");
        run("Select None");
        getStatistics(area, meanV, minV, maxV, std);
    } else {
    requireWindow(maskTitle);
    setOption("BlackBackground", true);
    setThreshold(1, 255);
    run("Create Selection");
    if (selectionType() == -1) {
        closeRegionsTempWindows();
        return 0;
    }
    roiManager("Add");
    activateImage("Regions_proc");
    roiManager("Select", 0);
    getStatistics(area, meanV, minV, maxV, std);
    roiManager("reset");
    run("Select None");
    }
    rangeOut[0] = minV;
    rangeOut[1] = maxV;
    rangeOut[2] = maxV - minV;
    if (rangeOut[2] <= 0) rangeOut[2] = 1;
    if (isFlatten) {
        if (peaksBright) regionResidualSpan = maxV - meanV;
        else regionResidualSpan = meanV - minV;
        if (regionResidualSpan < 0) regionResidualSpan = 0;
        // Bail only when there is genuinely no above/below-baseline excess at all.
        // (Do not bail on maxV<=meanV: a low-amplitude-but-valid residual is still usable.)
        if (peaksBright && maxV <= 0) {
            closeRegionsTempWindows();
            showStatus("No positive excess after flattening — try a larger flatten diameter or lower tolerance.");
            return 0;
        }
        if (!peaksBright && minV >= 0) {
            closeRegionsTempWindows();
            showStatus("No negative excess after flattening — try a larger flatten diameter or lower tolerance.");
            return 0;
        }
        thr = thresholdAbsFromMeanMax(meanV, minV, maxV, tolPct, peaksBright);
    } else {
        regionResidualSpan = maxV - minV;
        thr = thresholdAbsFromPct(minV, maxV, tolPct, peaksBright);
    }

    isIntensityWs = (indexOf(regionSepMode, "Intensity") >= 0);
    isBinaryWs = (indexOf(regionSepMode, "Binary") >= 0);

    // NOTE: do NOT multiply Regions_proc by the 0/255 mask here. That rescales the in-mask
    // intensities by 255 while thr/maxV were computed from the un-scaled image, so nothing
    // passes the threshold. Spatial restriction to the neurite is done by AND-ing the binary
    // result with the mask (below); the threshold runs on the true (un-scaled) residual.

    // Mask-restricted copy of the working image to seed the maxima watershed.
    // Clear outside the neurite mask (selection-based, no intensity rescaling) so the
    // watershed is seeded only by in-mask fluorescence and stray background maxima
    // do not create spurious ridge lines that cut into the patches.
    // In adjacent-particle mode the seed image is left un-cleared: off-mask patches must be
    // able to seed their own basins, which clearing outside the mask would prevent.
    if (isIntensityWs) {
        activateImage("Regions_proc");
        duplicateCurrentAs("Regions_intensity");
        if (!adjacentMode) {
        requireWindow(maskTitle);
        setOption("BlackBackground", true);
        setThreshold(1, 255);
        run("Create Selection");
        if (selectionType() != -1) {
            roiManager("reset");
            roiManager("Add");
            activateImage("Regions_intensity");
            roiManager("Select", 0);
            setBackgroundColor(0, 0, 0);
            run("Clear Outside");
            run("Select None");
            roiManager("reset");
        }
        activateImage(maskTitle);
        run("Select None");
        }
    }

    activateImage("Regions_proc");
    setOption("BlackBackground", true);
    if (peaksBright) setThreshold(thr, maxV);
    else setThreshold(minV, thr);
    run("Convert to Mask");
    rename("Regions_binary");
    // Flattening ran on the full image; restrict thresholded patches to the neurite mask.
    // Adjacent-particle mode keeps the full band: particles are rejected later by their
    // distance to the mask rather than clipped at its boundary.
    if (!adjacentMode) {
        imageCalculator("AND create", "Regions_binary", maskTitle);
        closeWindowIfOpen("Regions_binary");
        if (imageIsOpen("Result of Regions_binary")) { activateImage("Result of Regions_binary"); rename("Regions_binary"); }
    }

    // Remove sub-min-area noise speckles before watershed and counting.
    um2PerPxOpen = umPerPx * umPerPx;
    openPx = round(sqrt(minAreaUm2 / 3.14159265) / umPerPx);
    if (openPx < 1) openPx = 1;
    activateImage("Regions_binary");
    setOption("BlackBackground", true);
    for (oi = 0; oi < openPx; oi++)
        run("Open");

    // Save threshold-only mask for stage counting and fallback (all separation modes).
    separationWasApplied = 0;
    activateImage("Regions_binary");
    duplicateCurrentAs("Regions_binary_pre");
    regionCountThrFrag = countAllParticlesOnBinary("Regions_binary_pre", ch2Win);
    regionDiagN = captureAllRegionsFromBinary("Regions_binary_pre", ch2Win, regionDiagBXWork, regionDiagBYWork, regionDiagWWork, regionDiagHWork, regionDiagBXWork.length);
    regionCountThr = countParticlesOnBinary("Regions_binary_pre", ch2Win, minAreaUm2, maxAreaUm2, umPerPx);

    // ----- Separate touching/merged regions -----
    if (isIntensityWs) {
        separationWasApplied = 1;
        // Intensity-seeded watershed: each prominence-qualified local maximum becomes its own
        // basin; the watershed ridge lines cut the threshold mask so dense patches split apart.
        if (isFlatten) {
            if (peaksBright) rangeProc = maxV - meanV;
            else rangeProc = meanV - minV;
        } else {
            rangeProc = maxV - minV;
        }
        if (rangeProc <= 0) rangeProc = 1;
        promAbs = (regionPromPct / 100.0) * rangeProc;
        if (promAbs <= 0) promAbs = 1;
        fmOpt = "";
        if (!peaksBright) fmOpt = " light";
        activateImage("Regions_intensity");
        run("Find Maxima...", "prominence=" + d2s(promAbs, 6) + fmOpt + " output=[Segmented Particles]");
        segTitle = getTitle();
        imageCalculator("AND create", "Regions_binary", segTitle);
        closeWindowIfOpen("Regions_binary");
        if (imageIsOpen("Result of Regions_binary")) { activateImage("Result of Regions_binary"); rename("Regions_binary"); }
        closeWindowIfOpen(segTitle);
        closeWindowIfOpen("Regions_intensity");

        // Segmented Particles leaves label IDs 1..N; Analyze Particles needs 0/255 binary.
        binarizeNonzeroPixels("Regions_binary");
        activateImage("Regions_binary");
        getStatistics(areaWs, meanWs, minWs, maxWs, stdWs);
        if (maxWs <= 0) {
            closeWindowIfOpen("Regions_binary");
            activateImage("Regions_binary_pre");
            run("Select None");
            duplicateCurrentAs("Regions_binary");
            showStatus("No maxima at this prominence; using threshold mask only.");
        }
    } else if (isBinaryWs) {
        separationWasApplied = 1;
        // Classic distance-transform watershed on the binary region mask.
        activateImage("Regions_binary");
        setOption("BlackBackground", true);
        run("Watershed");
    }

    if (separationWasApplied) {
        regionCountFrag = countAllParticlesOnBinary("Regions_binary", ch2Win);
    } else {
        regionCountFrag = regionCountThr;
    }

    um2PerPx = umPerPx * umPerPx;
    minPx = round(minAreaUm2 / um2PerPx);
    if (minPx < 1) minPx = 1;
    if (maxAreaUm2 > 0) maxPx = round(maxAreaUm2 / um2PerPx);
    else maxPx = 99999999;
    sizeStr = "" + minPx + "-" + maxPx;
    measureTitle = "Regions_measure";
    prepareRegionsMeasureImage(ch2Win, measureTitle);
    // Intensity (Mean) must be read from the straightened analysis image, not the
    // 0/255 binary; the redirect target is only honored when set via Set Measurements.
    // Pin both the analyzed binary and the measure image to umPerPx so geometry is in microns.
    setMicronScale("Regions_binary", umPerPx);
    setMicronScale(measureTitle, umPerPx);
    run("Set Measurements...", "area mean bounding feret's shape redirect=" + measureTitle + " decimal=3");
    if (isOpen("Results")) { selectWindow("Results"); run("Clear Results"); }
    roiManager("reset");
    activateImage("Regions_binary");
    run("Analyze Particles...", "size=" + sizeStr + " pixel circularity=0.00-1.00 show=Nothing display clear add redirect=" + measureTitle);
    nOut = nResults;
    for (i = 0; i < nOut; i++) {
        peakRegBXWork[i] = getResult("BX", i);
        peakRegBYWork[i] = getResult("BY", i);
        peakRegWWork[i] = getResult("Width", i);
        peakRegHWork[i] = getResult("Height", i);
        xs[i] = peakRegBXWork[i] + peakRegWWork[i] / 2.0;
        ys[i] = peakRegBYWork[i] + peakRegHWork[i] / 2.0;
        peakRegMeanWork[i] = getResult("Mean", i);
        vals[i] = peakRegMeanWork[i];
        peakRegAreaWork[i] = getResult("Area", i);
        peakRegFeretWork[i] = getResult("Feret", i);
        peakRegMinFeretWork[i] = getResult("MinFeret", i);
        peakRegCircWork[i] = getResult("Circ.", i);
        peakRegSrcIdxWork[i] = i;
    }
    resetRegionMeasurementsRedirect();
    closeWindowIfOpen(measureTitle);
    if (isOpen("Results")) { selectWindow("Results"); run("Close"); }

    // Watershed can over-split into fragments below min area; revert to threshold mask.
    if (nOut == 0 && separationWasApplied) {
        if (regionCountThr > 0 || regionCountThrFrag > 0 || regionCountFrag > 0) {
            closeWindowIfOpen("Regions_binary");
            activateImage("Regions_binary_pre");
            run("Select None");
            duplicateCurrentAs("Regions_binary");
            if (isOpen("Results")) { selectWindow("Results"); run("Clear Results"); }
            roiManager("reset");
            prepareRegionsMeasureImage(ch2Win, measureTitle);
            setMicronScale("Regions_binary", umPerPx);
            setMicronScale(measureTitle, umPerPx);
            activateImage("Regions_binary");
            run("Set Measurements...", "area mean bounding feret's shape redirect=" + measureTitle + " decimal=3");
            activateImage("Regions_binary");
            run("Analyze Particles...", "size=" + sizeStr + " pixel circularity=0.00-1.00 show=Nothing display clear add redirect=" + measureTitle);
            nOut = nResults;
            for (i = 0; i < nOut; i++) {
                peakRegBXWork[i] = getResult("BX", i);
                peakRegBYWork[i] = getResult("BY", i);
                peakRegWWork[i] = getResult("Width", i);
                peakRegHWork[i] = getResult("Height", i);
                xs[i] = peakRegBXWork[i] + peakRegWWork[i] / 2.0;
                ys[i] = peakRegBYWork[i] + peakRegHWork[i] / 2.0;
                peakRegMeanWork[i] = getResult("Mean", i);
                vals[i] = peakRegMeanWork[i];
                peakRegAreaWork[i] = getResult("Area", i);
                peakRegFeretWork[i] = getResult("Feret", i);
                peakRegMinFeretWork[i] = getResult("MinFeret", i);
                peakRegCircWork[i] = getResult("Circ.", i);
                peakRegSrcIdxWork[i] = i;
            }
            resetRegionMeasurementsRedirect();
            closeWindowIfOpen(measureTitle);
            if (isOpen("Results")) { selectWindow("Results"); run("Close"); }
            if (nOut > 0) {
                showStatus("Watershed over-split all regions below min area; using threshold mask only.");
            }
        }
        if (nOut == 0 && regionCountThrFrag > 0) {
            showStatus(d2s(regionCountThrFrag, 0) + " speckles found but all below min area — lower Min region area or raise Tolerance.");
        }
    }

    // Measure each particle's distance to the neurite mask and midline, and drop those too far
    // away. Runs while the particle outlines are still in the ROI Manager (the previews and the
    // overlay export read them back through peakRegSrcIdxWork).
    regionCountPreDist = nOut;
    if (adjacentMode) nOut = applyAdjacentDistanceFilter(adjCoreMaskTitle, nOut, umPerPx, xs, ys, vals);

    regionCountWs = nOut;
    warnIfPatchCoordsNaN(xs, nOut);
    closeRegionsTempWindows();
    return nOut;
}

// Reduce the thresholded mask to the component(s) on the straightened midline (the traced
// neurite). Thresholding also picks up neighbouring structures elsewhere in the band; those
// fragments must not contribute to widths, in-mask metrics, particles, or adjacent distances.
// Leaves the result open as Mask_core (adjCoreMaskTitle) and fills maskFirstYPx / maskLastYPx /
// maskCorePresent from it.
function buildCoreMask() {
    adjCoreMaskTitle = maskTitle;
    requireWindow(maskTitle);
    run("Select None");
    getDimensions(cmW, cmH, cmC, cmS, cmF);
    midRow = round(maskMidYPx);
    if (midRow < 0) midRow = 0;
    if (midRow > cmH - 1) midRow = cmH - 1;
    // A few pixels of tolerance, so a slightly off-centre trace still finds its own neurite.
    seedHalfPx = round(cmH * 0.05);
    if (seedHalfPx < 1) seedHalfPx = 1;
    seedTop = midRow - seedHalfPx;
    if (seedTop < 0) seedTop = 0;
    seedBot = midRow + seedHalfPx;
    if (seedBot > cmH - 1) seedBot = cmH - 1;

    closeWindowIfOpen("Mask_seed");
    closeWindowIfOpen("Mask_core");
    activateImage(maskTitle);
    duplicateCurrentAs("Mask_seed");
    setColor(0);
    if (seedTop > 0) { makeRectangle(0, 0, cmW, seedTop); fill(); }
    if (seedBot < cmH - 1) { makeRectangle(0, seedBot + 1, cmW, cmH - seedBot - 1); fill(); }
    run("Select None");

    roiManager("reset");
    if (isOpen("Results")) { selectWindow("Results"); run("Clear Results"); }
    run("Set Measurements...", "area redirect=None decimal=3");
    activateImage(maskTitle);
    run("Analyze Particles...", "size=0-Infinity circularity=0.00-1.00 show=Nothing display clear add");
    nComp = nResults;
    compArea = newArray(nComp + 1);
    for (i = 0; i < nComp; i++) compArea[i] = getResult("Area", i);
    if (isOpen("Results")) { selectWindow("Results"); run("Close"); }
    resetRegionMeasurementsRedirect();

    if (nComp < 1) {
        roiManager("reset");
        closeWindowIfOpen("Mask_seed");
        showStatus("Core mask: thresholded mask is empty; keeping Mask_final.");
        return;
    }

    keepFlag = newArray(nComp);
    keepCount = 0;
    bestArea = -1;
    bestIdx = 0;
    for (i = 0; i < nComp; i++) {
        activateImage("Mask_seed");
        roiManager("Select", i);
        getStatistics(aSeed, mSeed, minSeed, maxSeed, sdSeed);
        run("Select None");
        if (maxSeed > 0) { keepFlag[i] = 1; keepCount++; } else keepFlag[i] = 0;
        if (compArea[i] > bestArea) { bestArea = compArea[i]; bestIdx = i; }
    }
    if (keepCount == 0) {
        // Nothing reaches the midline band (badly off-centre trace); fall back to the largest piece.
        keepFlag[bestIdx] = 1;
        keepCount = 1;
        showStatus("Core mask: no mask component on the midline; using the largest component.");
    }

    activateImage(maskTitle);
    run("Select None");
    duplicateCurrentAs("Mask_core");
    idMaskCore = pinCurrentImage();
    setColor(0);
    for (i = 0; i < nComp; i++) {
        if (keepFlag[i] == 0) {
            activateImage("Mask_core");
            roiManager("Select", i);
            fill();
        }
    }
    activateImage("Mask_core");
    run("Select None");
    roiManager("reset");
    closeWindowIfOpen("Mask_seed");
    adjCoreMaskTitle = "Mask_core";

    // Per column, keep only the contiguous run that contains (or lies nearest to) the midline, so a
    // branching or hooked component cannot contribute an edge far from the neurite.
    maskCorePresent = newArray(cmW);
    activateImage(adjCoreMaskTitle);
    for (x = 0; x < cmW; x++) {
        runStart = -1;
        if (getPixel(x, midRow) > 0) {
            runStart = midRow;
        } else {
            upY = -1;
            for (y = midRow - 1; y >= 0; y--) {
                if (getPixel(x, y) > 0) { upY = y; break; }
            }
            dnY = -1;
            for (y = midRow + 1; y < cmH; y++) {
                if (getPixel(x, y) > 0) { dnY = y; break; }
            }
            if (upY < 0) runStart = dnY;
            else if (dnY < 0) runStart = upY;
            else if ((midRow - upY) <= (dnY - midRow)) runStart = upY;
            else runStart = dnY;
        }
        if (runStart < 0) {
            maskFirstYPx[x] = -1;
            maskLastYPx[x] = -1;
            maskCorePresent[x] = 0;
        } else {
            first = runStart;
            last = runStart;
            while (first - 1 >= 0 && getPixel(x, first - 1) > 0) first--;
            while (last + 1 < cmH && getPixel(x, last + 1) > 0) last++;
            maskFirstYPx[x] = first;
            maskLastYPx[x] = last;
            maskCorePresent[x] = 1;
        }
    }

    // Carry the extent across columns where the neurite mask breaks, so the distance stays defined.
    // maskCorePresent is left untouched, so the width profile can still tell the two apart.
    fillFirst = -1;
    fillLast = -1;
    for (x = 0; x < cmW; x++) {
        if (maskFirstYPx[x] >= 0) { fillFirst = maskFirstYPx[x]; fillLast = maskLastYPx[x]; }
        else if (fillFirst >= 0) { maskFirstYPx[x] = fillFirst; maskLastYPx[x] = fillLast; }
    }
    fillFirst = -1;
    fillLast = -1;
    for (x = cmW - 1; x >= 0; x--) {
        if (maskFirstYPx[x] >= 0) { fillFirst = maskFirstYPx[x]; fillLast = maskLastYPx[x]; }
        else if (fillFirst >= 0) { maskFirstYPx[x] = fillFirst; maskLastYPx[x] = fillLast; }
    }
}

// Y-axis gap in pixels between the row span [pTop, pBot] and the mask foreground in column x.
// 0 when they overlap; -1 when the mask has no foreground in that column.
function maskGapPxInColumn(x, pTop, pBot) {
    if (x < 0 || x >= maskFirstYPx.length) return -1;
    mTop = maskFirstYPx[x];
    if (mTop < 0) return -1;
    mBot = maskLastYPx[x];
    if (pBot < mTop) return mTop - pBot;
    if (pTop > mBot) return pTop - mBot;
    return 0;
}

// Y-axis gap in pixels between the row span [pTop, pBot] and a single row; 0 when it falls inside.
function rowSpanGapPx(pTop, pBot, row) {
    if (row < pTop) return pTop - row;
    if (row > pBot) return row - pBot;
    return 0;
}

// Adjacent-particle mode: measure each particle's Y-axis distance to the neurite mask (nearest edge
// and centroid) and to the straightened-image midline, then reject those farther than
// adjMaxDistMidUm from the midline. All three distances are recorded either way.
// Whole particles are dropped rather than clipped, so area and diameter measurements stay intact.
// Compacts every peakReg*Work array in place and returns the surviving count; peakRegSrcIdxWork
// keeps each survivor pointing at its ROI Manager outline.
function applyAdjacentDistanceFilter(coreTitle, nOut, umPerPx, xs, ys, vals) {
    if (nOut < 1) return 0;
    requireWindow(coreTitle);
    nRoi = roiManager("count");
    kept = 0;
    midRowPx = round(maskMidYPx);
    for (i = 0; i < nOut; i++) {
        mgrIdx = peakRegSrcIdxWork[i];
        if (mgrIdx >= 0 && mgrIdx < nRoi) {
            activateImage(coreTitle);
            roiManager("Select", mgrIdx);
            getStatistics(aSel, mSel, minSel, maxSel, sdSel);
            getSelectionBounds(bx, by, bw, bh);
            run("Select None");
            overlapsMask = (maxSel > 0);
        } else {
            // No outline available: fall back to the calibrated bounding box.
            bx = round(peakRegBXWork[i] / umPerPx);
            by = round(peakRegBYWork[i] / umPerPx);
            bw = round(peakRegWWork[i] / umPerPx);
            bh = round(peakRegHWork[i] / umPerPx);
            overlapsMask = false;
        }
        if (bw < 1) bw = 1;
        if (bh < 1) bh = 1;
        pTop = by;
        pBot = by + bh - 1;
        cxPx = round(bx + bw / 2.0);
        cyPx = by + bh / 2.0;

        if (overlapsMask) {
            distMaskPx = 0;
        } else {
            distMaskPx = -1;
            for (x = bx; x < bx + bw; x++) {
                gap = maskGapPxInColumn(x, pTop, pBot);
                if (gap >= 0) {
                    if (distMaskPx < 0 || gap < distMaskPx) distMaskPx = gap;
                }
            }
        }
        distCentPx = maskGapPxInColumn(cxPx, round(cyPx), round(cyPx));
        // Gap filling in buildCoreMask leaves the mask extent defined in every column, so
        // this only fires when the mask holds no foreground at all; measure from the midline then.
        if (distMaskPx < 0) distMaskPx = rowSpanGapPx(pTop, pBot, midRowPx);
        if (distCentPx < 0) distCentPx = rowSpanGapPx(round(cyPx), round(cyPx), midRowPx);
        distMidUm = (cyPx - maskMidYPx) * umPerPx;

        keepParticle = (abs(distMidUm) <= adjMaxDistMidUm);

        if (keepParticle) {
            xs[kept] = xs[i];
            ys[kept] = ys[i];
            vals[kept] = vals[i];
            peakRegAreaWork[kept] = peakRegAreaWork[i];
            peakRegMeanWork[kept] = peakRegMeanWork[i];
            peakRegBXWork[kept] = peakRegBXWork[i];
            peakRegBYWork[kept] = peakRegBYWork[i];
            peakRegWWork[kept] = peakRegWWork[i];
            peakRegHWork[kept] = peakRegHWork[i];
            peakRegFeretWork[kept] = peakRegFeretWork[i];
            peakRegMinFeretWork[kept] = peakRegMinFeretWork[i];
            peakRegCircWork[kept] = peakRegCircWork[i];
            peakRegSrcIdxWork[kept] = peakRegSrcIdxWork[i];
            peakRegDistMaskWork[kept] = distMaskPx * umPerPx;
            peakRegDistMaskCentWork[kept] = distCentPx * umPerPx;
            peakRegDistMidWork[kept] = distMidUm;
            kept++;
        }
    }
    activateImage(coreTitle);
    run("Select None");
    return kept;
}

function warnIfPatchCoordsNaN(xs, nOut) {
    if (nOut < 1) return;
    if (!isNaN(xs[0])) return;
    print("Neurite FluoroQuant warning: patch bounding-box coords are NaN — Set Measurements must include bounding.");
    showStatus("Warning: patch positions invalid (bounding box not measured).");
}

// Like applyPeakSeparation but also filters peakReg*Work globals.
function applyPeakSeparationWithRegions(xs, ys, vals, nIn, minSepUm, umPerPx) {
    if (nIn < 1) return 0;
    // xs/ys come from Analyze Particles on a micron-calibrated image, so the separation has to be
    // compared in microns too. With no lower floor, 0 genuinely disables the filter.
    sepUm = minSepUm;
    if (sepUm < 0) sepUm = 0;
    idxOrder = sortIndicesByValueDesc(vals, nIn);
    keptX = newArray(nIn);
    keptY = newArray(nIn);
    keptV = newArray(nIn);
    keptRegArea = newArray(nIn);
    keptRegMean = newArray(nIn);
    keptRegBX = newArray(nIn);
    keptRegBY = newArray(nIn);
    keptRegW = newArray(nIn);
    keptRegH = newArray(nIn);
    keptRegFeret = newArray(nIn);
    keptRegMinFeret = newArray(nIn);
    keptRegCirc = newArray(nIn);
    keptRegSrcIdx = newArray(nIn);
    keptRegDistMask = newArray(nIn);
    keptRegDistMaskCent = newArray(nIn);
    keptRegDistMid = newArray(nIn);
    kCount = 0;
    for (ii = 0; ii < nIn; ii++) {
        idx = idxOrder[ii];
        xi = xs[idx];
        yi = ys[idx];
        vi = vals[idx];
        accept = 1;
        for (jj = 0; jj < kCount; jj++) {
            dxs = xi - keptX[jj];
            dys = yi - keptY[jj];
            if (sepUm > 0 && sqrt(dxs * dxs + dys * dys) < sepUm) {
                accept = 0;
                break;
            }
        }
        if (accept) {
            keptX[kCount] = xi;
            keptY[kCount] = yi;
            keptV[kCount] = vi;
            // peakRegSrcIdxWork is the identity unless adjacent-particle mode compacted the
            // arrays, so the outline in the ROI Manager is found through it rather than idx.
            peakRegionKeptRoiIdx[kCount] = peakRegSrcIdxWork[idx];
            keptRegSrcIdx[kCount] = peakRegSrcIdxWork[idx];
            keptRegDistMask[kCount] = peakRegDistMaskWork[idx];
            keptRegDistMaskCent[kCount] = peakRegDistMaskCentWork[idx];
            keptRegDistMid[kCount] = peakRegDistMidWork[idx];
            keptRegArea[kCount] = peakRegAreaWork[idx];
            keptRegMean[kCount] = peakRegMeanWork[idx];
            keptRegBX[kCount] = peakRegBXWork[idx];
            keptRegBY[kCount] = peakRegBYWork[idx];
            keptRegW[kCount] = peakRegWWork[idx];
            keptRegH[kCount] = peakRegHWork[idx];
            keptRegFeret[kCount] = peakRegFeretWork[idx];
            keptRegMinFeret[kCount] = peakRegMinFeretWork[idx];
            keptRegCirc[kCount] = peakRegCircWork[idx];
            kCount++;
        }
    }
    for (i = 0; i < kCount; i++) {
        xs[i] = keptX[i];
        ys[i] = keptY[i];
        vals[i] = keptV[i];
        peakRegAreaWork[i] = keptRegArea[i];
        peakRegMeanWork[i] = keptRegMean[i];
        peakRegBXWork[i] = keptRegBX[i];
        peakRegBYWork[i] = keptRegBY[i];
        peakRegWWork[i] = keptRegW[i];
        peakRegHWork[i] = keptRegH[i];
        peakRegFeretWork[i] = keptRegFeret[i];
        peakRegMinFeretWork[i] = keptRegMinFeret[i];
        peakRegCircWork[i] = keptRegCirc[i];
        peakRegSrcIdxWork[i] = keptRegSrcIdx[i];
        peakRegDistMaskWork[i] = keptRegDistMask[i];
        peakRegDistMaskCentWork[i] = keptRegDistMaskCent[i];
        peakRegDistMidWork[i] = keptRegDistMid[i];
    }
    return kCount;
}

function straightXPxToArcUm(xPx, arcWidthPx, roiLenUm) {
    if (arcWidthPx > 1) return (xPx / (arcWidthPx - 1)) * roiLenUm;
    return 0;
}

function sortIndicesByValueAsc(vals, n) {
    idx = newArray(n);
    for (i = 0; i < n; i++) idx[i] = i;
    for (i = 0; i < n - 1; i++) {
        best = i;
        for (j = i + 1; j < n; j++) {
            if (vals[idx[j]] < vals[idx[best]]) best = j;
        }
        tmp = idx[i];
        idx[i] = idx[best];
        idx[best] = tmp;
    }
    return idx;
}

function meanOfArray(vals, n) {
    if (n < 1) return 0;
    s = 0;
    for (i = 0; i < n; i++) s += vals[i];
    return s / n;
}

// Population CV (SD / mean) using only entries whose flag is non-zero. Empty or
// single-sample series, and non-positive means, return 0.
function cvOfValidValues(values, validFlags, n) {
    if (n < 1) return 0;
    sumN = 0;
    sumsqN = 0;
    nValid = 0;
    for (i = 0; i < n; i++) {
        if (validFlags[i]) {
            v = values[i];
            sumN += v;
            sumsqN += v * v;
            nValid++;
        }
    }
    if (nValid < 2) return 0;
    meanCv = sumN / nValid;
    if (meanCv <= 0) return 0;
    varCv = (sumsqN / nValid) - (meanCv * meanCv);
    if (varCv < 0) varCv = 0;
    return sqrt(varCv) / meanCv;
}

function meanAbsOfArray(vals, n) {
    if (n < 1) return 0;
    s = 0;
    for (i = 0; i < n; i++) s += abs(vals[i]);
    return s / n;
}

function meanPatchEquivDiamUm(regArea, n) {
    if (n < 1) return 0;
    s = 0;
    for (i = 0; i < n; i++) {
        if (regArea[i] > 0) s += 2.0 * sqrt(regArea[i] / 3.14159265);
    }
    return s / n;
}

function csvField(val) {
    if (val == "") return "";
    s = "" + val;
    if (indexOf(s, ",") >= 0 || indexOf(s, "\"") >= 0) {
        s = replace(s, "\"", "\"\"");
        return "\"" + s + "\"";
    }
    return s;
}

function masterSummaryRunId(base, outDirPath, autoIncrementOutput) {
    folderName = outDirPath;
    while (endsWith(folderName, "/") || endsWith(folderName, "\\")) {
        folderName = substring(folderName, 0, lengthOf(folderName) - 1);
    }
    if (autoIncrementOutput) return File.getName(folderName);
    return base + "_output";
}

// Slim master CSV: only columns that the current mode / enabled metrics actually fill.
// Use a separate master file per analysis mode (and preferably the same metric options).
function buildMasterSummaryHeader() {
    masterCsvHeaderStr = "File name";
    if (adjacentMode) {
        masterCsvHeaderStr = masterCsvHeaderStr + ",Analysis Ch1 Patch Count,Analysis Ch1 Patch Density per (um2),Analysis Ch1 Patch Mean Area (um2),Analysis Ch1 Patch Mean Intensity (AU),Analysis Ch1 Patch Maximum Diameter Mean (um),Analysis Ch1 Patch Minimum Diameter Mean (um),Analysis Ch1 Patch Circularity Mean,Analysis Ch1 Patch Equivalent Circle Diameter Mean (um)";
        masterCsvHeaderStr = masterCsvHeaderStr + ",Analysis Ch1 Patch Linear Density (patches / um),Analysis Ch1 Patch Mean Distance From Patch Edge To Mask (um),Analysis Ch1 Patch Mean Distance From Patch Centre To Midline (um)";
        if (adjDualMode) {
            masterCsvHeaderStr = masterCsvHeaderStr + ",Analysis Ch2 Patch Count,Analysis Ch2 Patch Density per (um2),Analysis Ch2 Patch Mean Area (um2),Analysis Ch2 Patch Mean Intensity (AU),Analysis Ch2 Patch Maximum Diameter Mean (um),Analysis Ch2 Patch Minimum Diameter Mean (um),Analysis Ch2 Patch Circularity Mean,Analysis Ch2 Patch Equivalent Circle Diameter Mean (um),Analysis Ch2 Patch Linear Density (patches / um),Analysis Ch2 Patch Mean Distance From Patch Edge To Mask (um),Analysis Ch2 Patch Mean Distance From Patch Centre To Midline (um)";
        }
    } else {
        if (includeCh2IntSummary) {
            masterCsvHeaderStr = masterCsvHeaderStr + ",Analysis Ch Mean Intensity Raw (AU),Analysis Ch Mean Intensity Bgsub (AU),Analysis Ch Projection Integrated Intensity Raw (AU),Analysis Ch Projection Integrated Intensity Bgsub (AU)";
        }
        if (includeCvSummary) {
            masterCsvHeaderStr = masterCsvHeaderStr + ",Analysis Channel Intensity CV";
            if (doDetrend) masterCsvHeaderStr = masterCsvHeaderStr + ",Analysis Channel Intensity CV Detrended";
        }
        if (doPeakDetection && isRegionMethodUsed) {
            masterCsvHeaderStr = masterCsvHeaderStr + ",Analysis Ch Patch Count,Analysis Ch Patch Density per (um2),Analysis Ch Patch Mean Area (um2),Analysis Ch Patch Mean Intensity (AU),Analysis Ch Patch Maximum Diameter Mean (um),Analysis Ch Patch Minimum Diameter Mean (um),Analysis Ch Patch Circularity Mean,Analysis Ch Patch Equivalent Circle Diameter Mean (um)";
        }
        if (doPeakDetection && !isRegionMethodUsed) {
            masterCsvHeaderStr = masterCsvHeaderStr + ",Analysis Ch Peak Count,Analysis Ch Peak Density per (um2)";
        }
        if (includeRoiSummary) {
            masterCsvHeaderStr = masterCsvHeaderStr + ",ROI Ch Length (um),ROI Ch Mean Intensity Raw (AU),ROI Ch Mean Intensity Bgsub (AU),ROI Ch Projection Integrated Intensity Raw (AU),ROI Ch Projection Integrated Intensity Bgsub (AU),ROI Channel Intensity CV,Neurite Width CV";
        }
        if (includeRoiSummary && includeMaskWidth) {
            masterCsvHeaderStr = masterCsvHeaderStr + ",Neurite Width Mean (um)";
        }
    }
    return masterCsvHeaderStr;
}

function buildMasterSummaryRow() {
    runId = masterSummaryRunId(base, outDir, autoIncrementOutput);
    line = csvField(runId);
    if (adjacentMode) {
        if (adjMasterReady) {
            line = line + "," + csvField(adjMasterNPeaks1) + "," + csvField(adjMasterPeakDensity1) + "," + csvField(adjMasterMeanArea1) + "," + csvField(adjMasterMeanInt1) + "," + csvField(adjMasterMeanFeret1) + "," + csvField(adjMasterMeanMinFeret1) + "," + csvField(adjMasterMeanCirc1) + "," + csvField(adjMasterMeanEquiv1);
            line = line + "," + csvField(adjMasterLinear1) + "," + csvField(adjMasterDistMask1) + "," + csvField(adjMasterDistMid1);
            if (adjDualMode) {
                line = line + "," + csvField(adjMasterNPeaks2) + "," + csvField(adjMasterPeakDensity2) + "," + csvField(adjMasterMeanArea2) + "," + csvField(adjMasterMeanInt2) + "," + csvField(adjMasterMeanFeret2) + "," + csvField(adjMasterMeanMinFeret2) + "," + csvField(adjMasterMeanCirc2) + "," + csvField(adjMasterMeanEquiv2) + "," + csvField(adjMasterLinear2) + "," + csvField(adjMasterDistMask2) + "," + csvField(adjMasterDistMid2);
            }
        } else {
            line = line + "," + csvField(nPeaks) + "," + csvField(peakDensity) + "," + csvField(meanOfArray(regAreaExport, nPeaks)) + "," + csvField(meanOfArray(regMeanExport, nPeaks)) + "," + csvField(meanOfArray(regFeretExport, nPeaks)) + "," + csvField(meanOfArray(regMinFeretExport, nPeaks)) + "," + csvField(meanOfArray(regCircExport, nPeaks)) + "," + csvField(meanPatchEquivDiamUm(regAreaExport, nPeaks));
            line = line + "," + csvField(adjLinearDensity) + "," + csvField(meanOfArray(regDistMaskExport, nPeaks)) + "," + csvField(meanAbsOfArray(regDistMidExport, nPeaks));
            if (adjDualMode) {
                line = line + ",,,,,,,,,,,";
            }
        }
    } else {
        if (includeCh2IntSummary) {
            line = line + "," + csvField(ch2_mean_inMask_raw) + "," + csvField(ch2_mean_inMask_bgsub) + "," + csvField(projIntC2Raw) + "," + csvField(projIntC2Bg);
        }
        if (includeCvSummary) {
            line = line + "," + csvField(cv_global);
            if (doDetrend) line = line + "," + csvField(cv_detrended);
        }
        if (doPeakDetection && isRegionMethodUsed) {
            line = line + "," + csvField(nPeaks) + "," + csvField(peakDensity) + "," + csvField(meanOfArray(regAreaExport, nPeaks)) + "," + csvField(meanOfArray(regMeanExport, nPeaks)) + "," + csvField(meanOfArray(regFeretExport, nPeaks)) + "," + csvField(meanOfArray(regMinFeretExport, nPeaks)) + "," + csvField(meanOfArray(regCircExport, nPeaks)) + "," + csvField(meanPatchEquivDiamUm(regAreaExport, nPeaks));
        }
        if (doPeakDetection && !isRegionMethodUsed) {
            line = line + "," + csvField(nPeaks) + "," + csvField(peakDensity);
        }
        if (includeRoiSummary) {
            line = line + "," + csvField(roiLenUm) + "," + csvField(meanC1Raw) + "," + csvField(meanC1Bg) + "," + csvField(projIntC1Raw) + "," + csvField(projIntC1Bg) + "," + csvField(roiCvGlobal) + "," + csvField(maskWidthCv);
        }
        if (includeRoiSummary && includeMaskWidth) {
            line = line + "," + csvField(maskWidthMeanUm);
        }
    }
    return line;
}

function masterSummaryFileIsBlank(path) {
    if (!File.exists(path)) return true;
    len = File.length(path);
    if (len == 0) return true;
    // Tiny leftover files (e.g. only a newline from an earlier blank-header write) count as blank.
    if (len < 16) {
        existing = File.openAsString(path);
        existing = replace(existing, "\r", "");
        existing = replace(existing, "\n", "");
        existing = replace(existing, "\t", "");
        existing = replace(existing, " ", "");
        if (existing == "") return true;
    }
    return false;
}

function appendMasterSummaryRow(path, dataLine) {
    if (path == "") return;
    parentDir = File.getDirectory(path);
    if (parentDir != "" && !File.exists(parentDir)) {
        showStatus("Master summary CSV: parent folder does not exist.");
        return;
    }
    // Blank / missing: overwrite with header + first data row in one write (ImageJ-safe).
    // Non-blank: append data only.
    if (masterSummaryFileIsBlank(path)) {
        masterCsvHeaderStr = buildMasterSummaryHeader();
        if (lengthOf(masterCsvHeaderStr) < 20) {
            showMessage("Master summary CSV", "Header string was empty or too short; not writing master CSV.");
            return;
        }
        File.saveString(masterCsvHeaderStr + "\n" + dataLine + "\n", path);
    } else {
        File.append(dataLine, path);
    }
}

function savePatchesTable(nPatches, regArea, regMean, regCirc, xs, ys, regBX, regBY, regW, regH, regFeret, regMinFeret, straightWidthPx, roiLenUm, umPerPx, outPath) {
    Table.create("Patches");
    if (nPatches < 1) {
        Table.save(outPath);
        if (isOpen("Patches")) { selectWindow("Patches"); run("Close"); }
        return;
    }
    idxOrder = sortIndicesByValueAsc(xs, nPatches);
    for (outRow = 0; outRow < nPatches; outRow++) {
        i = idxOrder[outRow];
        // Region geometry from Analyze Particles is already calibrated in microns (the straightened
        // image carries umPerPx). Report directly; re-applying a pixel->um conversion would double-scale.
        centreXUm = xs[i];
        centreYUm = ys[i];
        bboxMinXUm = regBX[i];
        bboxWidthXUm = regW[i];
        if (bboxWidthXUm < 0) bboxWidthXUm = 0;
        bboxMinYUm = regBY[i];
        bboxHeightYUm = regH[i];
        if (regArea[i] > 0) equivDiam = 2.0 * sqrt(regArea[i] / 3.14159265);
        else equivDiam = 0;
        Table.set("Patch number", outRow, outRow + 1);
        Table.set("Patch area (um2)", outRow, regArea[i]);
        Table.set("Mean fluorescence intensity (AU)", outRow, regMean[i]);
        Table.set("Maximum patch diameter (um)", outRow, regFeret[i]);
        Table.set("Minimum patch diameter (um)", outRow, regMinFeret[i]);
        Table.set("Circularity", outRow, regCirc[i]);
        Table.set("Equivalent circle diameter (um)", outRow, equivDiam);
        Table.set("Patch centre X (um along neurite)", outRow, centreXUm);
        Table.set("Patch centre Y (um across neurite)", outRow, centreYUm);
        Table.set("BBox min X (um along neurite)", outRow, bboxMinXUm);
        Table.set("BBox min Y (um across neurite)", outRow, bboxMinYUm);
        Table.set("BBox width X (um along neurite)", outRow, bboxWidthXUm);
        Table.set("BBox height Y (um across neurite)", outRow, bboxHeightYUm);
        if (adjacentMode) {
            // Headings must not contain commas: Table.save writes them unquoted, so a comma would
            // split one heading into two CSV fields and misalign the whole header row.
            Table.set("Distance from patch edge to mask (um)", outRow, regDistMaskExport[i]);
            Table.set("Distance from patch centre to mask (um)", outRow, regDistMaskCentExport[i]);
            Table.set("Distance from patch centre to midline (um)", outRow, abs(regDistMidExport[i]));
        }
    }
    Table.save(outPath);
    if (isOpen("Patches")) { selectWindow("Patches"); run("Close"); }
}

// Indices sorted by vals descending (selection sort; avoids float tie issues in rankPositions).
function sortIndicesByValueDesc(vals, n) {
    idx = newArray(n);
    for (i = 0; i < n; i++) idx[i] = i;
    for (i = 0; i < n - 1; i++) {
        best = i;
        for (j = i + 1; j < n; j++) {
            if (vals[idx[j]] > vals[idx[best]]) best = j;
        }
        tmp = idx[i];
        idx[i] = idx[best];
        idx[best] = tmp;
    }
    return idx;
}

// Sort values ascending (selection sort by value).
function sortValuesAscending(src, n) {
    sorted = newArray(n);
    if (n < 1) return sorted;
    idx = newArray(n);
    for (i = 0; i < n; i++) idx[i] = i;
    for (i = 0; i < n - 1; i++) {
        best = i;
        for (j = i + 1; j < n; j++) {
            if (src[idx[j]] < src[idx[best]]) best = j;
        }
        tmp = idx[i];
        idx[i] = idx[best];
        idx[best] = tmp;
    }
    for (i = 0; i < n; i++) sorted[i] = src[idx[i]];
    return sorted;
}

function mapArcIndicesToPeaks(peakIdx, profileArc, yAtMax, xs, ys, vals) {
    nIdx = peakIdx.length;
    for (i = 0; i < nIdx; i++) {
        px = peakIdx[i];
        xs[i] = px;
        ys[i] = yAtMax[px];
        vals[i] = profileArc[px];
    }
    return nIdx;
}

function applyPeakSeparation(xs, ys, vals, nIn, minSepUm, umPerPx) {
    if (nIn < 1) return 0;
    sepPx = minSepUm / umPerPx;
    if (sepPx < 1) sepPx = 1;
    sepPx = round(sepPx);
    idxOrder = sortIndicesByValueDesc(vals, nIn);
    keptX = newArray(nIn);
    keptY = newArray(nIn);
    keptV = newArray(nIn);
    kCount = 0;
    for (ii = 0; ii < nIn; ii++) {
        idx = idxOrder[ii];
        xi = xs[idx];
        yi = ys[idx];
        vi = vals[idx];
        accept = 1;
        for (jj = 0; jj < kCount; jj++) {
            if (abs(xi - keptX[jj]) < sepPx) {
                accept = 0;
                break;
            }
        }
        if (accept) {
            keptX[kCount] = xi;
            keptY[kCount] = yi;
            keptV[kCount] = vi;
            kCount++;
        }
    }
    for (i = 0; i < kCount; i++) {
        xs[i] = keptX[i];
        ys[i] = keptY[i];
        vals[i] = keptV[i];
    }
    return kCount;
}

// peakCounts[0]=after FM/particles, [1]=after separation, [2]=final; fills xs, ys, vals (sized nIn then nFinal).
// Masked regions method also fills peakReg*Work globals.
function runPeakDetection(method, maskTitle, ch2Win, swArc, shArc, tolPct, minSepUm, peaksBright, doRefine, umPerPx, regionMinAreaUm, regionMaxAreaUm, regionSmoothUm, profileArc, yAtMax, arcValid, xs, ys, vals, peakCounts) {
    peakIdx = newArray(0);
    nFm = 0;
    isRegions = 0;
    if (indexOf(method, "Masked regions") >= 0) isRegions = 1;
    if (isRegions) {
        nFm = findRegions2DMasked(maskTitle, ch2Win, tolPct, peaksBright, regionMinAreaUm, regionMaxAreaUm, regionSmoothUm, umPerPx, xs, ys, vals);
    } else {
        tolAbsArc = computeArcTolAbs(profileArc, arcValid, swArc, tolPct);
        peakIdx = findPeaksArcFiji(profileArc, arcValid, swArc, tolAbsArc, peaksBright);
        nFm = peakIdx.length;
        if (nFm > 0) mapArcIndicesToPeaks(peakIdx, profileArc, yAtMax, xs, ys, vals);
    }
    peakCounts[0] = nFm;
    if (nFm < 1) {
        peakCounts[1] = 0;
        peakCounts[2] = 0;
        return;
    }
    if (isRegions) {
        nSep = applyPeakSeparationWithRegions(xs, ys, vals, nFm, minSepUm, umPerPx);
    } else {
        if (doRefine) {
            refinePx = 1;
            refinePeaksInMask(xs, ys, vals, nFm, maskTitle, ch2Win, swArc, shArc, refinePx);
        }
        nSep = applyPeakSeparation(xs, ys, vals, nFm, minSepUm, umPerPx);
    }
    peakCounts[1] = nSep;
    peakCounts[2] = nSep;
}

// Open peak preview at default (100%) zoom; user can adjust with + / - while dialog is open.
function zoomPeakPreviewForTuning(previewTitle) {
    peakPreviewZoomPct = 100;
    hasWin = 0;
    titles = getList("image.titles");
    for (i = 0; i < titles.length; i++) {
        if (startsWith(titles[i], previewTitle)) hasWin = 1;
    }
    if (hasWin == 0) return;
    selectWindowByTitlePrefix(previewTitle);
    run("Select None");
    run("Set... ", "zoom=100");
}

function buildPeakPreview(ch2Win, xs, ys, nPeaks, markerUm, umPerPx, previewTitle) {
    closeWindowIfOpen(previewTitle);
    selectWindow(ch2Win);
    run("Select None");
    run("Duplicate...", "title=" + previewTitle);
    selectWindow(previewTitle);
    run("RGB Color");
    run("Select None");
    setForegroundColor(255, 255, 0);
    setLineWidth(2);
    radPx = round((markerUm / 2.0) / umPerPx);
    if (radPx < 1) radPx = 1;
    for (i = 0; i < nPeaks; i++) {
        makeOval(xs[i] - radPx, ys[i] - radPx, 2 * radPx, 2 * radPx);
        run("Draw");
        run("Select None");
    }
    clearImageRoi();
}

function buildRegionPreview(ch2Win, nRegions, xs, ys, regBX, regBY, regW, regH, previewTitle) {
    closeWindowIfOpen(previewTitle);
    selectWindow(ch2Win);
    run("Select None");
    run("Duplicate...", "title=" + previewTitle);
    selectWindow(previewTitle);
    run("RGB Color");
    run("Select None");
    drawMaskOutlineOnWindow(previewTitle);
    drawAdjacentGuideLinesOnWindow(previewTitle);
    selectWindow(previewTitle);
    setForegroundColor(255, 255, 0);
    setLineWidth(2);
    nRoi = roiManager("count");
    drewRoi = 0;
    if (nRoi > 0) {
        for (i = 0; i < nRegions; i++) {
            mgrRoiIdx = peakRegionKeptRoiIdx[i];
            if (mgrRoiIdx >= 0) {
                if (mgrRoiIdx < nRoi) {
                    selectWindow(previewTitle);
                    roiManager("Select", mgrRoiIdx);
                    run("Draw");
                    run("Select None");
                    drewRoi = 1;
                }
            }
        }
    }
    if (drewRoi == 0) {
        for (i = 0; i < nRegions; i++) {
            makeRectangle(regBX[i], regBY[i], regW[i], regH[i]);
            run("Draw");
            run("Select None");
        }
    }
    if (drewRoi == 1 && nRegions > 0) {
        saveKeptRegionExportRois(nRegions);
    } else if (drewRoi == 0 && nRegions > 0) {
        selectWindow(ch2Win);
        roiManager("reset");
        for (i = 0; i < nRegions; i++) {
            makeRectangle(regBX[i], regBY[i], regW[i], regH[i]);
            roiManager("Add");
            run("Select None");
            peakRegionKeptRoiIdx[i] = i;
        }
        saveKeptRegionExportRois(nRegions);
    }
    roiManager("reset");
    clearImageRoi();
}

// 1 if this pixel bounding-box centre is inside the adjacent midline band, or if the band
// does not apply (non-adjacent, or scale not set). Same test as applyAdjacentDistanceFilter.
function speckleCentreWithinMidlineBand(by, bh) {
    if (!adjacentMode) return 1;
    if (adjUmPerPx <= 0) return 1;
    cyPx = by + bh / 2.0;
    distMidUm = (cyPx - maskMidYPx) * adjUmPerPx;
    if (abs(distMidUm) <= adjMaxDistMidUm) return 1;
    return 0;
}

// Orange outlines for sub-min-area speckles when no regions pass the area filter.
// Adjacent mode hides speckles whose centre lies outside the green midline band.
function buildRegionDiagnosticPreview(ch2Win, nRegions, regBX, regBY, regW, regH, previewTitle) {
    closeWindowIfOpen(previewTitle);
    selectWindow(ch2Win);
    run("Select None");
    run("Duplicate...", "title=" + previewTitle);
    selectWindow(previewTitle);
    run("RGB Color");
    run("Select None");
    drawMaskOutlineOnWindow(previewTitle);
    drawAdjacentGuideLinesOnWindow(previewTitle);
    selectWindow(previewTitle);
    setForegroundColor(255, 165, 0);
    setLineWidth(1);
    loadedDiagRois = 0;
    // Prefer the true pixel outlines saved during detection (irregular shapes), so the
    // diagnostic preview matches the actual speckles rather than bounding boxes.
    roiManager("reset");
    if (regionDiagRoiPath != "" && File.exists(regionDiagRoiPath)) {
        roiManager("Open", regionDiagRoiPath);
        nRoi = roiManager("count");
        if (nRoi > 0) loadedDiagRois = 1;
        for (i = 0; i < nRoi; i++) {
            selectWindow(previewTitle);
            roiManager("Select", i);
            getSelectionBounds(dbx, dby, dbw, dbh);
            if (speckleCentreWithinMidlineBand(dby, dbh) == 0) {
                run("Select None");
                continue;
            }
            run("Draw");
            run("Select None");
        }
        roiManager("reset");
    }
    if (loadedDiagRois == 0) {
        // Fallback: regBX/regBY/regW/regH come from getSelectionBounds (always pixel units),
        // and makeRectangle also expects pixels, so use them directly with no scaling.
        for (i = 0; i < nRegions; i++) {
            if (speckleCentreWithinMidlineBand(regBY[i], regH[i]) == 0) continue;
            makeRectangle(regBX[i], regBY[i], regW[i], regH[i]);
            run("Draw");
            run("Select None");
        }
    }
    roiManager("reset");
    clearImageRoi();
}

function buildArcProfilePlot(profileArc, arcValid, swArc, xs, nPeaks, plotTitle, roiLenUm) {
    closeWindowIfOpen(plotTitle);
    plotH = 80;
    marginL = 82;
    marginB = 54;
    if (plotH < 32) plotH = 32;
    if (swArc < 1) return;
    if (plotH < 1) return;
    plotW = swArc + marginL;
    plotTotalH = plotH + marginB;
    newImage(plotTitle, "8-bit", plotW, plotTotalH, 1);
    selectWindow(plotTitle);
    minP = 0;
    maxP = 1;
    hasP = 0;
    for (x = 0; x < swArc; x++) {
        if (arcValid[x]) {
            if (hasP == 0) {
                minP = profileArc[x];
                maxP = profileArc[x];
                hasP = 1;
            } else {
                if (profileArc[x] < minP) minP = profileArc[x];
                if (profileArc[x] > maxP) maxP = profileArc[x];
            }
        }
    }
    if (hasP == 0) return;
    rangeP = maxP - minP;
    if (rangeP <= 0) rangeP = 1;
    baseY = plotH - 1;
    ox = marginL;
    for (x = 0; x < plotW; x++) {
        for (y = 0; y < plotTotalH; y++) setPixel(x, y, 0);
    }
    prevPy = -1;
    for (x = 0; x < swArc; x++) {
        if (arcValid[x]) {
            normV = round((profileArc[x] - minP) / rangeP * 255);
            if (normV < 0) normV = 0;
            if (normV > 255) normV = 255;
            py = baseY - round((profileArc[x] - minP) / rangeP * (plotH - 1));
            if (py < 0) py = 0;
            if (py > baseY) py = baseY;
            pxPlot = ox + x;
            for (y = py; y <= baseY; y++) setPixel(pxPlot, y, normV);
            if (prevPy >= 0) {
                y0 = prevPy;
                y1 = py;
                if (y1 < y0) {
                    yt = y0;
                    y0 = y1;
                    y1 = yt;
                }
                for (y = y0; y <= y1; y++) setPixel(pxPlot, y, normV);
            }
            prevPy = py;
        } else {
            prevPy = -1;
        }
    }
    peakMark = 255;
    for (i = 0; i < nPeaks; i++) {
        px = round(xs[i]);
        if (px >= 0) {
            if (px < swArc) {
                if (arcValid[px]) {
                    py = baseY - round((profileArc[px] - minP) / rangeP * (plotH - 1));
                    if (py < 0) py = 0;
                    if (py > baseY) py = baseY;
                    pxPlot = ox + px;
                    for (dy = -3; dy <= 3; dy++) {
                        xx = pxPlot + dy;
                        if (xx >= ox) {
                            if (xx < ox + swArc) setPixel(xx, py, peakMark);
                        }
                    }
                    for (dy = -2; dy <= 2; dy++) {
                        yy = py + dy;
                        if (yy >= 0) {
                            if (yy < plotH) setPixel(pxPlot, yy, peakMark);
                        }
                    }
                }
            }
        }
    }
    setFont("SansSerif", 11, "antialiased");
    setColor(235, 235, 235);
    drawString("max", 6, 10);
    drawString(d2s(maxP, 2), 6, 22);
    xLblY = plotH + 12;
    drawString("min " + d2s(minP, 2), 6, xLblY);
    drawString("0", ox + 4, xLblY);
    xEndLbl = ox + swArc - 52;
    if (xEndLbl > plotW - 58) xEndLbl = plotW - 58;
    if (xEndLbl < ox + 4) xEndLbl = ox + 4;
    drawString(d2s(roiLenUm, 2) + " um", xEndLbl, xLblY);
    xMidLbl = ox + round(swArc / 2) - 48;
    if (xMidLbl < ox + 4) xMidLbl = ox + 4;
    drawString("Arc length (um)", xMidLbl, plotH + 26);
    setFont("SansSerif", 9, "antialiased");
    drawString("Analysis bgsub max in mask per column", 6, plotH + 38);
    drawString("Y: normalized display (raw values in CSV)", 6, plotTotalH - 8);
    resetMinAndMax();
}

function layoutPeakPreviews() {
    titles = getList("image.titles");
    hasPeakPrev = 0;
    for (i = 0; i < titles.length; i++) {
        if (startsWith(titles[i], peakPreviewPrefix)) hasPeakPrev = 1;
    }
    if (hasPeakPrev == 0) return;
    padX = 40;
    padY = 80;
    selectWindowByTitlePrefix(peakPreviewPrefix);
    winW = round(getWidth() * peakPreviewZoomPct / 100) + 28;
    winH = round(getHeight() * peakPreviewZoomPct / 100) + 28;
    if (winH < 220) winH = 220;
    setLocation(padX, padY, winW, winH);
    pw = winW;
    if (isOpen(peakProfilePlotWin)) {
        selectWindow(peakProfilePlotWin);
        setLocation(padX + pw + 20, padY);
        selectWindowByTitlePrefix(peakPreviewPrefix);
        clearImageRoi();
    }
}

// 3D tracing overlay: remove guidance marker before/after per-plane tracing.
function clearPathTraceOverlay(srcWin) {
    if (!isOpen(srcWin)) return;
    selectWindow(srcWin);
    Overlay.remove();
    run("Select None");
}

// 3D tracing overlay: green endpoint marker visible on all Z planes (overlay position 0).
function showLastEndpointMarker(srcWin, endX, endY) {
    if (!isOpen(srcWin)) return;
    selectWindow(srcWin);
    Overlay.remove();
    cx = round(endX);
    cy = round(endY);
    r = 6;
    setColor(0, 255, 0);
    setLineWidth(1);
    Overlay.drawLine(cx - r, cy - r, cx + r, cy + r);
    Overlay.setPosition(0);
    Overlay.drawLine(cx - r, cy + r, cx + r, cy - r);
    Overlay.setPosition(0);
    Overlay.show();
    run("Select None");
}

// Mean intensity of entire strip image (for per-segment Z compensation).
function meanStripIntensity(stripTitle) {
    selectWindow(stripTitle);
    getStatistics(area, mean, min, max, std);
    return mean;
}

// Scale strip intensities in place.
function scaleStripIntensity(stripTitle, factor) {
    if (factor <= 0 || (factor > 0.9999 && factor < 1.0001)) return;
    selectWindow(stripTitle);
    run("Multiply...", "value=" + factor);
}

// Snap segment start to previous segment end before straightening (reduces junction gap).
function snapSegmentStartToPreviousEnd(segIdx, prevEndX, prevEndY) {
    roiManager("Select", segIdx);
    rtSnap = selectionType();
    if (rtSnap != 5 && rtSnap != 6) return;
    getSelectionCoordinates(pxSnap, pySnap);
    if (pxSnap.length < 2) return;
    pxSnap[0] = prevEndX;
    pySnap[0] = prevEndY;
    if (rtSnap == 5) makeSelection("polyline", pxSnap, pySnap);
    else makeSelection("freeline", pxSnap, pySnap);
    roiManager("Update");
}

// Load runtime settings from def_* globals (batch mode skips setup dialogs).
function applySettingsFromDefs() {
    roiIdx = round(def_roiIdx);
    ch2Idx = round(def_ch2Idx);
    fovUm = def_fovUm;
    straightenUm = def_straightUm;
    tissueType = def_tissueType;
    compensateZIntensityPref = def_compensateZIntensity;
    is3DProject = (indexOf(tissueType, "Z-range") >= 0);
    is3DPerPlane = (indexOf(tissueType, "per-plane") >= 0);
    is3D = is3DProject || is3DPerPlane;
    if (is3DPerPlane) compensateZIntensity = compensateZIntensityPref;
    else compensateZIntensity = false;
    maskIdx = roiIdx;
    zProjMethod = def_zProjMethod;
    rbUm = def_rbUm;
    maskThr = def_maskThr;
    maskDilUm = def_maskDilUm;
    includeCh2IntSummary = def_includeCh2IntSummary;
    includeCh2Profile = def_includeCh2Profile;
    includeNP = def_includeNP;
    includeCvSummary = def_includeCvSummary;
    includeRoiSummary = def_includeRoiSummary;
    includeRoiProfile = def_includeRoiProfile;
    includeMaskWidth = def_includeMaskWidth;
    includeAnalysisParams = def_includeAnalysisParams;
    includeCv = includeCvSummary;
    wcvUm = def_wcvUm;
    doDetrend = def_doDetrend;
    detrendUm = def_detrendUm;
    includeCh1 = includeRoiSummary || includeRoiProfile;
    adjacentMode = def_adjacentMode;
    adjSecondAnalysisIdx = round(def_adjSecondAnalysisIdx);
    if (adjSecondAnalysisIdx < 0) adjSecondAnalysisIdx = 0;
    adjDualMode = 0;
    if (adjacentMode && adjSecondAnalysisIdx > 0 && adjSecondAnalysisIdx != ch2Idx) adjDualMode = 1;
    else adjSecondAnalysisIdx = 0;
    adjMaxDistMidUm = def_adjMaxDistMidUm;
    if (adjMaxDistMidUm < 0) adjMaxDistMidUm = 0;
}

// Pick Z range once (from ROI channel) and projection method; sets globals projZFirst, projZLast, zProjMethod.
function pickZRangeGlobal(srcWin) {
    requireWindow(srcWin);
    selectWindow(srcWin);
    if (!batchMode) {
        Dialog.create("Z-range projection");
        Dialog.addChoice("Projection method:", newArray("Max Intensity", "Average Intensity", "Sum Slices"), zProjMethod);
        Dialog.show();
        zProjMethod = Dialog.getChoice();
    } else {
        zProjMethod = def_zProjMethod;
    }
    waitForUser("Z range",
        "Scroll to the FIRST slice containing the neurite, then click OK.");
    zFirstPick = getSliceNumber();
    waitForUser("Z range",
        "Scroll to the LAST slice containing the neurite, then click OK.");
    zLastPick = getSliceNumber();
    if (zFirstPick > zLastPick) {
        tmpZ = zFirstPick;
        zFirstPick = zLastPick;
        zLastPick = tmpZ;
    }
    projZFirst = zFirstPick;
    projZLast = zLastPick;
    prefSet("unified.projZFirst", projZFirst);
    prefSet("unified.projZLast", projZLast);
}

// Z-project slices zFirst..zLast of srcWin and replace that window with the 2D projection.
function projectChannelStack(srcWin, zFirst, zLast, projMethod) {
    requireWindow(srcWin);
    if (zFirst > zLast) {
        tmpZ = zFirst;
        zFirst = zLast;
        zLast = tmpZ;
    }
    tempDup = srcWin + "_zsub";
    closeWindowIfOpen(tempDup);
    selectWindow(srcWin);
    if (nSlices > 1) run("Duplicate...", "title=" + tempDup + " duplicate");
    else run("Duplicate...", "title=" + tempDup);
    selectWindow(tempDup);
    run("Make Substack...", "channels=1 slices=" + zFirst + "-" + zLast + " frames=1");
    subTitle = getTitle();
    if (projMethod == "Average Intensity") projArg = "[Average Intensity]";
    else if (projMethod == "Sum Slices") projArg = "[Sum Slices]";
    else projArg = "[Max Intensity]";
    run("Z Project...", "projection=" + projArg);
    projName = getTitle();
    closeWindowIfOpen(srcWin);
    if (isOpen(tempDup) && tempDup != projName) closeWindowIfOpen(tempDup);
    if (isOpen(subTitle) && subTitle != projName && subTitle != tempDup) closeWindowIfOpen(subTitle);
    selectWindow(projName);
    rename(srcWin);
    run("Select None");
}

// 3D: non-blocking tracing dialog; user draws while dialog is open. Returns 1 to finish with
// existing segments (checkbox, nSeg>0 only), 0 to accept a new segment.
function waitForSegmentOrFinish(srcWin, nSeg) {
    selectWindow(srcWin);
    Dialog.createNonBlocking("3D path tracing");
    if (nSeg > 0) {
        Dialog.addMessage("Segment " + (nSeg + 1) + ": draw a connecting segmented line ROI on the next z-plane.\n" +
            "Green marker = end of previous segment.");
        Dialog.addCheckbox("Finish path (use segments already traced)", false);
    } else {
        Dialog.addMessage("Segment 1: draw a segmented line ROI on the first z-plane of interest.");
    }
    Dialog.show();
    if (nSeg > 0 && Dialog.getCheckbox()) return 1;
    return 0;
}

// 3D: capture the neurite as an ordered series of per-plane segments. Each segment is a
// Segmented/Freehand line traced on whichever slice it is in focus; the same plane may be
// reused any number of times. Fills the ROI Manager (one ROI per segment, in trace order),
// the global segZ[] (slice per segment) and pathTotalPx (summed in-plane arc length in px).
function capturePathSegments(srcWin) {
    roiManager("Reset");
    segZ = newArray(0);
    segStartX = newArray(0);
    segStartY = newArray(0);
    segEndX = newArray(0);
    segEndY = newArray(0);
    pathTotalPx = 0;
    nSeg = 0;
    clearPathTraceOverlay(srcWin);
    moreSegments = true;
    while (moreSegments) {
        selectWindow(srcWin);
        run("Select None");
        if (nSeg > 0) {
            showLastEndpointMarker(srcWin, segEndX[nSeg - 1], segEndY[nSeg - 1]);
        }
        finishRequested = waitForSegmentOrFinish(srcWin, nSeg);
        if (finishRequested == 1) {
            moreSegments = false;
            break;
        }
        rtSeg = selectionType();
        if (rtSeg != 5 && rtSeg != 6) {
            if (nSeg > 0) {
                showMessage("3D path tracing",
                    "No valid line drawn.\nUncheck \"Finish path\" and draw a segment, or check it and click OK.");
                continue;
            }
            cancelRunAndExit(activeOutDir, "");
        }
        getSelectionCoordinates(sx, sy);
        if (sx.length < 2) cancelRunAndExit(activeOutDir, "Each segment must contain at least two points.");
        zSeg = getSliceNumber();
        roiManager("Add");
        segZ = Array.concat(segZ, zSeg);
        segStartX = Array.concat(segStartX, sx[0]);
        segStartY = Array.concat(segStartY, sy[0]);
        segEndX = Array.concat(segEndX, sx[sx.length - 1]);
        segEndY = Array.concat(segEndY, sy[sy.length - 1]);
        for (i = 1; i < sx.length; i++) {
            dxS = sx[i] - sx[i-1]; dyS = sy[i] - sy[i-1];
            pathTotalPx += sqrt(dxS*dxS + dyS*dyS);
        }
        nSeg++;
        showLastEndpointMarker(srcWin, segEndX[nSeg - 1], segEndY[nSeg - 1]);
        Dialog.create("Add another segment?");
        Dialog.addMessage(nSeg + " segment(s) traced so far.");
        Dialog.addChoice("Action:", newArray("Add another plane segment", "Done - finish path"), "Add another plane segment");
        Dialog.show();
        actionSeg = Dialog.getChoice();
        moreSegments = (indexOf(actionSeg, "Add") >= 0);
    }
    clearPathTraceOverlay(srcWin);
    pathSegCount = nSeg;
    return nSeg;
}

// 3D: straighten each captured segment on its own slice and concatenate the strips left-to-right
// into a single 2D window (outWin), so downstream analysis sees the same kind of straightened image
// as the 2D path. Uses the segment ROIs already in the ROI Manager and the global segZ[].
function straightenAlongPath(srcWin, outWin, straightWidthPx, nSeg) {
    stripTitles = newArray(nSeg);
    totalW = 0;
    maxH = 0;
    bd = 8;
    for (i = 0; i < nSeg; i++) {
        selectWindow(srcWin);
        if (i > 0) snapSegmentStartToPreviousEnd(i, segEndX[i - 1], segEndY[i - 1]);
        roiManager("Select", i);
        Stack.setSlice(segZ[i]);
        run("Straighten...", "line=" + straightWidthPx);
        stripT = "PathStrip_" + outWin + "_" + i;
        rename(stripT);
        stripTitles[i] = stripT;
        totalW += getWidth();
        if (getHeight() > maxH) maxH = getHeight();
        bd = bitDepth();
    }
    if (maxH < 1) maxH = straightWidthPx;
    if (compensateZIntensity && nSeg > 1) {
        meanRef = 0;
        refIdx = 0;
        stripMeans = newArray(nSeg);
        for (i = 0; i < nSeg; i++) {
            stripMeans[i] = meanStripIntensity(stripTitles[i]);
            if (stripMeans[i] > meanRef) {
                meanRef = stripMeans[i];
                refIdx = i;
            }
        }
        for (i = 0; i < nSeg; i++) {
            if (i == refIdx) continue;
            meanI = stripMeans[i];
            if (meanI > 0 && meanRef > 0) scaleStripIntensity(stripTitles[i], meanRef / meanI);
        }
    }
    if (bd == 24) typeStr = "RGB black";
    else if (bd == 16) typeStr = "16-bit black";
    else if (bd == 32) typeStr = "32-bit black";
    else typeStr = "8-bit black";
    if (isOpen(outWin)) { selectWindow(outWin); run("Close"); }
    newImage(outWin, typeStr, totalW, maxH, 1);
    xOff = 0;
    for (i = 0; i < nSeg; i++) {
        selectWindow(stripTitles[i]);
        wStrip = getWidth();
        hStrip = getHeight();
        run("Select All");
        run("Copy");
        selectWindow(outWin);
        makeRectangle(xOff, 0, wStrip, hStrip);
        run("Paste");
        run("Select None");
        xOff += wStrip;
        selectWindow(stripTitles[i]);
        run("Close");
    }
    selectWindow(outWin);
    run("Select None");
    resetMinAndMax();
}

// Load detection parameters for the current adjacent pass (primary vs adj2 prefs).
function loadAdjPassPrefs(useAdj2) {
    if (useAdj2) {
        tolPct = def_adj2_tolPct;
        minSepUm = def_adj2_minSepUm;
        regionMinAreaUm = def_adj2_regionMinAreaUm;
        regionMaxAreaUm = def_adj2_regionMaxAreaUm;
        regionSmoothUm = def_adj2_regionSmoothUm;
        regionPromPct = def_adj2_regionPromPct;
    } else {
        tolPct = def_tolPct;
        minSepUm = def_minSepUm;
        regionMinAreaUm = def_regionMinAreaUm;
        regionMaxAreaUm = def_regionMaxAreaUm;
        regionSmoothUm = def_regionSmoothUm;
        regionPromPct = def_regionPromPct;
    }
    // Adjacent mode always uses absolute intensity (no local-baseline flatten).
    regionDogUm = 0;
    adjMaxDistMidUm = def_adjMaxDistMidUm;
    if (adjMaxDistMidUm < 0) adjMaxDistMidUm = 0;
    markerUm = def_markerUm;
    doRefine = def_doRefine;
}

function saveAdjPassPrefs(useAdj2) {
    if (useAdj2) {
        prefSet("unified.adj2.tolPct", tolPct);
        prefSet("unified.adj2.minSepUm", minSepUm);
        prefSet("unified.adj2.regionMinAreaUm", regionMinAreaUm);
        prefSet("unified.adj2.regionMaxAreaUm", regionMaxAreaUm);
        prefSet("unified.adj2.regionSmoothUm", regionSmoothUm);
        prefSet("unified.adj2.regionPromPct", regionPromPct);
        def_adj2_tolPct = tolPct;
        def_adj2_minSepUm = minSepUm;
        def_adj2_regionMinAreaUm = regionMinAreaUm;
        def_adj2_regionMaxAreaUm = regionMaxAreaUm;
        def_adj2_regionSmoothUm = regionSmoothUm;
        def_adj2_regionPromPct = regionPromPct;
    } else {
        prefSet("unified.tolPct", tolPct);
        prefSet("unified.minSepUm", minSepUm);
        prefSet("unified.regionMinAreaUm", regionMinAreaUm);
        prefSet("unified.regionMaxAreaUm", regionMaxAreaUm);
        prefSet("unified.regionSmoothUm", regionSmoothUm);
        prefSet("unified.regionPromPct", regionPromPct);
        def_tolPct = tolPct;
        def_minSepUm = minSepUm;
        def_regionMinAreaUm = regionMinAreaUm;
        def_regionMaxAreaUm = regionMaxAreaUm;
        def_regionSmoothUm = regionSmoothUm;
        def_regionPromPct = regionPromPct;
    }
    prefSet("unified.adjacentMaxDistMidUm", adjMaxDistMidUm);
    def_adjMaxDistMidUm = adjMaxDistMidUm;
}

// Point analysisCorrWin / export tags / preview titles at the channel for this adjacent pass.
function setAdjacentPassContext(chIdx, corrWin, exportTag, winBase) {
    analysisCorrWin = corrWin;
    analysisExportTag = exportTag;
    masterAnalysisChIdx = chIdx;
    analysisOverlayWin = winBase + "_overlay";
    peakPreviewPrefix = winBase + "_peak_preview";
    peakProfilePlotWin = winBase + "_peak_profile";
    analysisPlainWin = winBase + "_plain";
}

// Snapshot patch metrics and detection tuning for the master / parameters CSVs.
function stashAdjacentMasterMetrics() {
    meanArea = meanOfArray(regAreaExport, nPeaks);
    meanInt = meanOfArray(regMeanExport, nPeaks);
    meanFeret = meanOfArray(regFeretExport, nPeaks);
    meanMinFeret = meanOfArray(regMinFeretExport, nPeaks);
    meanCirc = meanOfArray(regCircExport, nPeaks);
    meanEquiv = meanPatchEquivDiamUm(regAreaExport, nPeaks);
    meanDistMask = meanOfArray(regDistMaskExport, nPeaks);
    meanDistMid = meanAbsOfArray(regDistMidExport, nPeaks);
    if (adjPassUseAdj2Prefs) {
        adjMasterChIdx2 = masterAnalysisChIdx;
        adjMasterNPeaks2 = nPeaks;
        adjMasterPeakDensity2 = peakDensity;
        adjMasterMeanArea2 = meanArea;
        adjMasterMeanInt2 = meanInt;
        adjMasterMeanFeret2 = meanFeret;
        adjMasterMeanMinFeret2 = meanMinFeret;
        adjMasterMeanCirc2 = meanCirc;
        adjMasterMeanEquiv2 = meanEquiv;
        adjMasterLinear2 = adjLinearDensity;
        adjMasterDistMask2 = meanDistMask;
        adjMasterDistMid2 = meanDistMid;
        adjParamTol2 = tolPct;
        adjParamMinSep2 = minSepUm;
        adjParamMinArea2 = regionMinAreaUm;
        adjParamMaxArea2 = regionMaxAreaUm;
        adjParamSmooth2 = regionSmoothUm;
        adjParamProm2 = regionPromPct;
        adjParamDog2 = regionDogUm;
    } else {
        adjMasterChIdx1 = masterAnalysisChIdx;
        adjMasterNPeaks1 = nPeaks;
        adjMasterPeakDensity1 = peakDensity;
        adjMasterMeanArea1 = meanArea;
        adjMasterMeanInt1 = meanInt;
        adjMasterMeanFeret1 = meanFeret;
        adjMasterMeanMinFeret1 = meanMinFeret;
        adjMasterMeanCirc1 = meanCirc;
        adjMasterMeanEquiv1 = meanEquiv;
        adjMasterLinear1 = adjLinearDensity;
        adjMasterDistMask1 = meanDistMask;
        adjMasterDistMid1 = meanDistMid;
        adjParamTol1 = tolPct;
        adjParamMinSep1 = minSepUm;
        adjParamMinArea1 = regionMinAreaUm;
        adjParamMaxArea1 = regionMaxAreaUm;
        adjParamSmooth1 = regionSmoothUm;
        adjParamProm1 = regionPromPct;
        adjParamDog1 = regionDogUm;
    }
}

// Patches, summary and overlay for the current adjacent analysis channel (master CSV is appended once later).
function exportAdjacentChannelResults() {
    totalMaskArea = 0;
    for (x = 0; x < widthsUm.length; x++) totalMaskArea += widthsUm[x] * umPerPx;
    if (totalMaskArea > 0) peakDensity = nPeaks / totalMaskArea; else peakDensity = 0;
    if (roiLenUm > 0) adjLinearDensity = nPeaks / roiLenUm; else adjLinearDensity = 0;

    savePatchesTable(nPeaks, regAreaExport, regMeanExport, regCircExport, xs, ys, regBXExport, regBYExport, regWExport, regHExport, regFeretExport, regMinFeretExport, swPk, roiLenUm, umPerPx, analysisOutDir + File.separator + base + analysisExportTag + "_patches.csv");

    Table.create("Analysis Summary");
    if (includeMaskWidth) Table.set("Neurite Width Mean (um)", 0, maskWidthMeanUm);
    Table.set("Analysis Patch Count", 0, nPeaks);
    Table.set("Analysis Patch Density per (um2)", 0, peakDensity);
    Table.set("Analysis Patch Mean Area (um2)", 0, meanOfArray(regAreaExport, nPeaks));
    Table.set("Analysis Patch Mean Intensity (AU)", 0, meanOfArray(regMeanExport, nPeaks));
    Table.set("Analysis Patch Maximum Diameter Mean (um)", 0, meanOfArray(regFeretExport, nPeaks));
    Table.set("Analysis Patch Minimum Diameter Mean (um)", 0, meanOfArray(regMinFeretExport, nPeaks));
    Table.set("Analysis Patch Circularity Mean", 0, meanOfArray(regCircExport, nPeaks));
    Table.set("Analysis Patch Equivalent Circle Diameter Mean (um)", 0, meanPatchEquivDiamUm(regAreaExport, nPeaks));
    Table.set("Analysis Patch Linear Density (patches / um)", 0, adjLinearDensity);
    Table.set("Analysis Patch Mean Distance From Patch Edge To Mask (um)", 0, meanOfArray(regDistMaskExport, nPeaks));
    Table.set("Analysis Patch Mean Distance From Patch Centre To Midline (um)", 0, meanAbsOfArray(regDistMidExport, nPeaks));
    Table.save(analysisOutDir + File.separator + base + analysisExportTag + "_summary.csv");
    if (isOpen("Analysis Summary")) { selectWindow("Analysis Summary"); run("Close"); }

    activateImage(analysisCorrWin);
    duplicateCurrentAs(analysisOverlayWin);
    run("Remove Overlay"); setOption("Display Label", false);
    addPeakOverlaysToImage(1, nPeaks, xs, ys, markerUm, umPerPx);
    run("Flatten");
    overlayPeaksPath = analysisOutDir + File.separator + base + analysisExportTag + "_overlay_patches.png";
    saveAs("PNG", overlayPeaksPath);
    closeWindowAfterSave(analysisOverlayWin, overlayPeaksPath);
    // Keep a durable copy of this pass's outlines for the dual merged overlay (the shared temp zip
    // is overwritten by the next pass).
    if (adjDualMode && regionExportRoiPath != "" && File.exists(regionExportRoiPath)) {
        zipDest = getDirectory("temp") + "neurite_adjRoi_C" + d2s(masterAnalysisChIdx, 0) + ".zip";
        if (File.exists(zipDest)) File.delete(zipDest);
        File.copy(regionExportRoiPath, zipDest);
        if (adjPassUseAdj2Prefs) adjRoiZipPass2 = zipDest;
        else adjRoiZipPass1 = zipDest;
    }
    roiManager("reset");
    stashAdjacentMasterMetrics();
}

// Interactive adjacent-particle tuning for the channel currently in analysisCorrWin.
function runAdjacentInteractiveTuningPass() {
    peakMethod = "Masked regions (threshold + particles)";
    requireWindow(maskTitle); getDimensions(swMpk, shMpk, cMpk, sMpk, fMpk);
    activateImage(analysisCorrWin);
    getDimensions(swPk, shPk, cPk, sPk, fPk);
    swArc = swPk;
    if (swMpk < swArc) swArc = swMpk;
    shArc = shPk;
    if (shMpk < shArc) shArc = shMpk;
    profileArc = newArray(swArc);
    yAtMax = newArray(swArc);
    arcValid = newArray(swArc);
    peakBufSize = swArc * shArc;
    if (peakBufSize < 1) peakBufSize = 1;
    xsWork = newArray(peakBufSize);
    ysWork = newArray(peakBufSize);
    valsWork = newArray(peakBufSize);
    peakRegAreaWork = newArray(peakBufSize);
    peakRegMeanWork = newArray(peakBufSize);
    peakRegBXWork = newArray(peakBufSize);
    peakRegBYWork = newArray(peakBufSize);
    peakRegWWork = newArray(peakBufSize);
    peakRegHWork = newArray(peakBufSize);
    peakRegFeretWork = newArray(peakBufSize);
    peakRegMinFeretWork = newArray(peakBufSize);
    peakRegCircWork = newArray(peakBufSize);
    peakRegSrcIdxWork = newArray(peakBufSize);
    peakRegDistMaskWork = newArray(peakBufSize);
    peakRegDistMaskCentWork = newArray(peakBufSize);
    peakRegDistMidWork = newArray(peakBufSize);
    peakCounts = newArray(3);
    peakRegionKeptRoiIdx = newArray(peakBufSize);
    regionDiagBXWork = newArray(peakBufSize);
    regionDiagBYWork = newArray(peakBufSize);
    regionDiagWWork = newArray(peakBufSize);
    regionDiagHWork = newArray(peakBufSize);
    peaksBright = true;
    regionExportRoiPath = "";
    setBatchMode(false);
    peakTuneTitle = "Patch detection tuning — Adjacent patches (C" + d2s(masterAnalysisChIdx, 0) + ")";
    donePeaks = false;
    peakPreviewReady = 0;
    while (!donePeaks) {
        Dialog.createNonBlocking(peakTuneTitle);
        Dialog.addNumber("Tolerance (% of range):", tolPct, 2, 5, "Min patch brightness; lower = more patches");
        Dialog.addNumber("Min inter-peak distance (um):", minSepUm, 2, 5, "Drop peaks closer than this");
        Dialog.addNumber("Min region area (um^2):", regionMinAreaUm, 3, 5, "Min patch size; lower = smaller blobs");
        Dialog.addNumber("Max region area (um^2, 0=none):", regionMaxAreaUm, 1, 5, "0 = no maximum patch size");
        Dialog.addNumber("Pre-threshold Gaussian sigma (um, 0=off):", regionSmoothUm, 2, 5, "Blur analysis channel before detect; lower = sharper");
        Dialog.addNumber("Maxima prominence (% of band range):", regionPromPct, 1, 5, "Split merged patches; lower = more splits");
        bandHalfUm = (shPk / 2.0) * umPerPx;
        if (adjMaxDistMidUm > bandHalfUm) distHint = "Exceeds the straightened band (" + d2s(bandHalfUm, 2) + " um each side)";
        else distHint = "Reject particles farther than this from the midline";
        Dialog.addNumber("Max distance from midline (um):", adjMaxDistMidUm, 2, 5, distHint);
        if (peakPreviewReady == 1) {
            if (nPeaks > 0) {
                totalRegArea = sumRegionAreas(nPeaks);
                meanRegArea = totalRegArea / nPeaks;
                countMsg = "Blobs: " + d2s(nPeaks, 0) +
                    "  |  Mean area: " + d2s(meanRegArea, 2) + " um^2" +
                    "  |  Total area: " + d2s(totalRegArea, 2) + " um^2";
                countMsg = countMsg + "\nMean distance to mask: " + d2s(meanOfArray(peakRegDistMaskWork, nPeaks), 2) +
                    " um  |  Mean |distance| to midline: " + d2s(meanAbsOfArray(peakRegDistMidWork, nPeaks), 2) + " um";
            } else {
                if (regionCountPreDist > 0) {
                    countMsg = "Blobs: 0 (" + d2s(regionCountPreDist, 0) + " passed size but all outside max distance from midline)";
                } else if (regionDiagN > 0) {
                    countMsg = "Blobs: 0 (orange = below min area inside the band)";
                } else {
                    countMsg = "Blobs: 0";
                }
            }
            Dialog.addMessage(countMsg);
        }
        if (peakPreviewReady == 1) peakActionDefault = "Accept patches and continue";
        else peakActionDefault = "Update preview";
        Dialog.addChoice("Action:", newArray("Update preview", "Accept patches and continue"), peakActionDefault);
        Dialog.show();
        tolPct = Dialog.getNumber();
        minSepUm = Dialog.getNumber();
        regionMinAreaUm = Dialog.getNumber();
        regionMaxAreaUm = Dialog.getNumber();
        regionSmoothUm = Dialog.getNumber();
        regionPromPct = Dialog.getNumber();
        adjMaxDistMidUm = Dialog.getNumber();
        if (adjMaxDistMidUm < 0) adjMaxDistMidUm = 0;
        regionDogUm = 0;
        actionPeaks = Dialog.getChoice();
        saveAdjPassPrefs(adjPassUseAdj2Prefs);
        doPreviewUpdate = 0;
        doAcceptPeaks = 0;
        if (indexOf(actionPeaks, "Update") >= 0) doPreviewUpdate = 1;
        if (indexOf(actionPeaks, "Accept") >= 0) {
            doAcceptPeaks = 1;
            if (peakPreviewReady == 0) {
                showStatus("Running detection before accept...");
                doPreviewUpdate = 1;
            }
        }
        if (doPreviewUpdate) {
            requireWindow(maskTitle);
            closeWindowsWithTitlePrefix(peakPreviewPrefix);
            closeWindowIfOpen(peakProfilePlotWin);
            roiManager("reset");
            beginHiddenCompute();
            runPeakDetection(peakMethod, maskTitle, analysisCorrWin, swArc, shArc, tolPct, minSepUm, peaksBright, doRefine, umPerPx, regionMinAreaUm, regionMaxAreaUm, regionSmoothUm, profileArc, yAtMax, arcValid, xsWork, ysWork, valsWork, peakCounts);
            endHiddenCompute();
            nPeaksArcFM = peakCounts[0];
            nPeaksArcSep = peakCounts[1];
            nPeaks = peakCounts[2];
            xs = newArray(nPeaks);
            ys = newArray(nPeaks);
            vals = newArray(nPeaks);
            for (i = 0; i < nPeaks; i++) {
                xs[i] = xsWork[i];
                ys[i] = ysWork[i];
                vals[i] = valsWork[i];
            }
            if (nPeaks > 0) {
                buildRegionPreview(analysisCorrWin, nPeaks, xs, ys, peakRegBXWork, peakRegBYWork, peakRegWWork, peakRegHWork, peakPreviewPrefix);
            } else if (regionCountPreDist > 0) {
                buildRegionPreview(analysisCorrWin, 0, xs, ys, peakRegBXWork, peakRegBYWork, peakRegWWork, peakRegHWork, peakPreviewPrefix);
            } else if (regionDiagN > 0) {
                buildRegionDiagnosticPreview(analysisCorrWin, regionDiagN, regionDiagBXWork, regionDiagBYWork, regionDiagWWork, regionDiagHWork, peakPreviewPrefix);
            } else {
                buildRegionPreview(analysisCorrWin, 0, xs, ys, peakRegBXWork, peakRegBYWork, peakRegWWork, peakRegHWork, peakPreviewPrefix);
            }
            zoomPeakPreviewForTuning(peakPreviewPrefix);
            layoutPeakPreviews();
            selectWindowByTitlePrefix(peakPreviewPrefix);
            peakPreviewReady = 1;
            if (nPeaks == 0 && regionCountPreDist > 0) {
                showStatus(d2s(regionCountPreDist, 0) + " particles passed size but all were outside the max distance from midline.");
            } else if (nPeaks == 0 && regionDiagN > 0) {
                showStatus("Diagnostic preview (orange): " + d2s(regionDiagN, 0) + " speckles below min area inside the band — raise Tolerance or lower Min region area.");
            }
        }
        if (doAcceptPeaks) {
            closeWindowsWithTitlePrefix(peakPreviewPrefix);
            closeWindowIfOpen(peakProfilePlotWin);
            isRegionMethodUsed = 1;
            detectionType = "regions";
            regAreaExport = newArray(nPeaks);
            regMeanExport = newArray(nPeaks);
            regBXExport = newArray(nPeaks);
            regBYExport = newArray(nPeaks);
            regWExport = newArray(nPeaks);
            regHExport = newArray(nPeaks);
            regFeretExport = newArray(nPeaks);
            regMinFeretExport = newArray(nPeaks);
            regCircExport = newArray(nPeaks);
            regDistMaskExport = newArray(nPeaks);
            regDistMaskCentExport = newArray(nPeaks);
            regDistMidExport = newArray(nPeaks);
            for (i = 0; i < nPeaks; i++) {
                regAreaExport[i] = peakRegAreaWork[i];
                regMeanExport[i] = peakRegMeanWork[i];
                regBXExport[i] = peakRegBXWork[i];
                regBYExport[i] = peakRegBYWork[i];
                regWExport[i] = peakRegWWork[i];
                regHExport[i] = peakRegHWork[i];
                regFeretExport[i] = peakRegFeretWork[i];
                regMinFeretExport[i] = peakRegMinFeretWork[i];
                regCircExport[i] = peakRegCircWork[i];
                regDistMaskExport[i] = peakRegDistMaskWork[i];
                regDistMaskCentExport[i] = peakRegDistMaskCentWork[i];
                regDistMidExport[i] = peakRegDistMidWork[i];
            }
            donePeaks = true;
        }
    }
}

// Dual adjacent only: RGB merge of ROI (R) + primary analysis (G) + second analysis (B), with
// primary patches in yellow and second-channel patches in cyan.
function buildAdjacentMergedOverlay() {
    if (!adjDualMode) return;
    if (adjMergeRoiWin == "" || !imageIsOpen(adjMergeRoiWin)) {
        showStatus("Merged overlay skipped: ROI merge source missing.");
        return;
    }
    if (!imageIsOpen(analysisCorrWin) || !imageIsOpen(analysis2CorrWin)) {
        showStatus("Merged overlay skipped: analysis channels missing.");
        return;
    }
    closeWindowIfOpen(adjMergedOverlayWin);
    closeWindowIfOpen("AdjMerge_R");
    closeWindowIfOpen("AdjMerge_G");
    closeWindowIfOpen("AdjMerge_B");
    activateImage(adjMergeRoiWin);
    duplicateCurrentAs("AdjMerge_R");
    resetMinAndMax();
    run("Enhance Contrast", "saturated=0.35");
    activateImage(analysisCorrWin);
    duplicateCurrentAs("AdjMerge_G");
    resetMinAndMax();
    run("Enhance Contrast", "saturated=0.35");
    activateImage(analysis2CorrWin);
    duplicateCurrentAs("AdjMerge_B");
    resetMinAndMax();
    run("Enhance Contrast", "saturated=0.35");
    run("Merge Channels...", "c1=AdjMerge_R c2=AdjMerge_G c3=AdjMerge_B create keep");
    mergeTitle = getTitle();
    if (mergeTitle != adjMergedOverlayWin) rename(adjMergedOverlayWin);
    activateImage(adjMergedOverlayWin);
    run("RGB Color");
    // RGB Color can rename the window on some builds.
    if (indexOf(getTitle(), adjMergedOverlayWin) != 0) rename(adjMergedOverlayWin);
    closeWindowIfOpen("AdjMerge_R");
    closeWindowIfOpen("AdjMerge_G");
    closeWindowIfOpen("AdjMerge_B");
    // Drop leftover composite titles Merge Channels may leave behind.
    closeWindowIfOpen("Composite");
    activateImage(adjMergedOverlayWin);
    run("Remove Overlay");
    setOption("Display Label", false);
    addMaskOutlineOverlay(adjMergedOverlayWin);
    activateImage(adjMergedOverlayWin);
    if (adjRoiZipPass1 != "" && File.exists(adjRoiZipPass1)) {
        roiManager("reset");
        roiManager("Open", adjRoiZipPass1);
        nExp1 = roiManager("count");
        for (i = 0; i < nExp1; i++) {
            roiManager("Select", i);
            run("Add Selection...", "stroke=yellow width=1");
            run("Select None");
        }
        roiManager("reset");
    }
    if (adjRoiZipPass2 != "" && File.exists(adjRoiZipPass2)) {
        roiManager("reset");
        roiManager("Open", adjRoiZipPass2);
        nExp2 = roiManager("count");
        for (i = 0; i < nExp2; i++) {
            roiManager("Select", i);
            run("Add Selection...", "stroke=cyan width=1");
            run("Select None");
        }
        roiManager("reset");
    }
    activateImage(adjMergedOverlayWin);
    if (Overlay.size > 0) Overlay.drawLabels(false);
    run("Flatten");
    mergedPath = analysisOutDir + File.separator + base + "_overlay_patches_merged.png";
    saveAs("PNG", mergedPath);
    setOption("Changes", false);
    closeWindowAfterSave(adjMergedOverlayWin, mergedPath);
    closeWindowsWithTitlePrefix(adjMergedOverlayWin);
}

// Adjacent mode only: one or two detection passes (primary, then optional second analysis channel).
function runAdjacentDetectionPasses() {
    doPeakDetection = 1;
    peakMethod = "Masked regions (threshold + particles)";
    isRegionMethodUsed = 1;
    detectionType = "regions";
    adjacentExportsDone = 0;
    adjMasterReady = 0;
    adjMasterChIdx2 = 0;
    adjMasterNPeaks2 = 0;
    adjMasterPeakDensity2 = 0;
    adjMasterMeanArea2 = 0;
    adjMasterMeanInt2 = 0;
    adjMasterMeanFeret2 = 0;
    adjMasterMeanMinFeret2 = 0;
    adjMasterMeanCirc2 = 0;
    adjMasterMeanEquiv2 = 0;
    adjMasterLinear2 = 0;
    adjMasterDistMask2 = 0;
    adjMasterDistMid2 = 0;
    adjParamTol1 = 0;
    adjParamMinSep1 = 0;
    adjParamMinArea1 = 0;
    adjParamMaxArea1 = 0;
    adjParamSmooth1 = 0;
    adjParamProm1 = 0;
    adjParamDog1 = 0;
    adjParamTol2 = 0;
    adjParamMinSep2 = 0;
    adjParamMinArea2 = 0;
    adjParamMaxArea2 = 0;
    adjParamSmooth2 = 0;
    adjParamProm2 = 0;
    adjParamDog2 = 0;
    totalMaskArea = 0;
    for (x = 0; x < widthsUm.length; x++) totalMaskArea += widthsUm[x] * umPerPx;

    adjPassCount = 1;
    if (adjDualMode) adjPassCount = 2;
    primaryCorrWin = analysisCorrWin;
    primaryExportTag = analysisExportTag;
    primaryWinBase = "Analysis_C" + d2s(ch2Idx, 0);

    for (adjPassIdx = 0; adjPassIdx < adjPassCount; adjPassIdx++) {
        adjPassUseAdj2Prefs = 0;
        if (adjPassIdx == 0) {
            setAdjacentPassContext(ch2Idx, primaryCorrWin, primaryExportTag, primaryWinBase);
        } else {
            adjPassUseAdj2Prefs = 1;
            setAdjacentPassContext(adjSecondAnalysisIdx, analysis2CorrWin, analysis2ExportTag, analysis2Win);
        }
        loadAdjPassPrefs(adjPassUseAdj2Prefs);
        regionExportRoiPath = "";
        if (batchMode) {
            runAutoPeakDetection();
        } else {
            runAdjacentInteractiveTuningPass();
        }
        requireWindow(maskTitle);
        beginHiddenCompute();
        finalizePeakDetectionResults();
        exportAdjacentChannelResults();
        closeWindowIfOpen(analysisOverlayWin);
        endHiddenCompute();
        if (nPeaks == 0 && !batchMode) {
            showMessage("Particle detection", "No particles found on C" + d2s(masterAnalysisChIdx, 0) + " at the current settings.\n" +
                "Try a lower tolerance, a smaller min area, or a larger max distance from midline.\n" +
                "After particles: " + d2s(nPeaksArcFM, 0) + " region(s).");
        }
    }
    // Restore primary context for any later shared references.
    setAdjacentPassContext(ch2Idx, primaryCorrWin, primaryExportTag, primaryWinBase);
    analysisCorrWin = primaryCorrWin;
    if (adjDualMode) {
        beginHiddenCompute();
        buildAdjacentMergedOverlay();
        closeWindowIfOpen(adjMergedOverlayWin);
        closeWindowIfOpen("AdjMerge_R");
        closeWindowIfOpen("AdjMerge_G");
        closeWindowIfOpen("AdjMerge_B");
        closeWindowIfOpen("Composite");
        endHiddenCompute();
    }
    // One master row per run: primary metrics in the usual columns, second channel in the trailing block.
    adjMasterReady = 1;
    masterAnalysisChIdx = 0;
    nPeaks = adjMasterNPeaks1;
    peakDensity = adjMasterPeakDensity1;
    adjLinearDensity = adjMasterLinear1;
    if (appendMasterSummary && masterSummaryPath != "") {
        masterRow = buildMasterSummaryRow();
        appendMasterSummaryRow(masterSummaryPath, masterRow);
    }
    adjMasterReady = 0;
    adjacentExportsDone = 1;
}

// Run peak detection using saved prefs (batch mode; no tuning dialogs).
function runAutoPeakDetection() {
    peakMethod = "Masked regions (threshold + particles)";
    prefSet("unified.peakMethod", peakMethod);
    if (adjacentMode && adjPassUseAdj2Prefs) {
        loadAdjPassPrefs(1);
    } else if (adjacentMode) {
        loadAdjPassPrefs(0);
    } else {
        tolPct = def_tolPct;
        minSepUm = def_minSepUm;
        markerUm = def_markerUm;
        doRefine = def_doRefine;
        regionMinAreaUm = def_regionMinAreaUm;
        regionMaxAreaUm = def_regionMaxAreaUm;
        regionSmoothUm = def_regionSmoothUm;
        regionPromPct = def_regionPromPct;
        regionDogUm = def_regionDogUm;
    }
    requireWindow(maskTitle);
    getDimensions(swMpk, shMpk, cMpk, sMpk, fMpk);
    activateImage(analysisCorrWin);
    getDimensions(swPk, shPk, cPk, sPk, fPk);
    swArc = swPk;
    if (swMpk < swArc) swArc = swMpk;
    shArc = shPk;
    if (shMpk < shArc) shArc = shMpk;
    profileArc = newArray(swArc);
    yAtMax = newArray(swArc);
    arcValid = newArray(swArc);
    peakBufSize = swArc * shArc;
    if (peakBufSize < 1) peakBufSize = 1;
    xsWork = newArray(peakBufSize);
    ysWork = newArray(peakBufSize);
    valsWork = newArray(peakBufSize);
    peakRegAreaWork = newArray(peakBufSize);
    peakRegMeanWork = newArray(peakBufSize);
    peakRegBXWork = newArray(peakBufSize);
    peakRegBYWork = newArray(peakBufSize);
    peakRegWWork = newArray(peakBufSize);
    peakRegHWork = newArray(peakBufSize);
    peakRegFeretWork = newArray(peakBufSize);
    peakRegMinFeretWork = newArray(peakBufSize);
    peakRegCircWork = newArray(peakBufSize);
    peakRegSrcIdxWork = newArray(peakBufSize);
    peakRegDistMaskWork = newArray(peakBufSize);
    peakRegDistMaskCentWork = newArray(peakBufSize);
    peakRegDistMidWork = newArray(peakBufSize);
    peakCounts = newArray(3);
    peakRegionKeptRoiIdx = newArray(peakBufSize);
    regionDiagBXWork = newArray(peakBufSize);
    regionDiagBYWork = newArray(peakBufSize);
    regionDiagWWork = newArray(peakBufSize);
    regionDiagHWork = newArray(peakBufSize);
    peaksBright = true;
    requireWindow(maskTitle);
    roiManager("reset");
    beginHiddenCompute();
    runPeakDetection(peakMethod, maskTitle, analysisCorrWin, swArc, shArc, tolPct, minSepUm, peaksBright, doRefine, umPerPx, regionMinAreaUm, regionMaxAreaUm, regionSmoothUm, profileArc, yAtMax, arcValid, xsWork, ysWork, valsWork, peakCounts);
    endHiddenCompute();
    nPeaksArcFM = peakCounts[0];
    nPeaksArcSep = peakCounts[1];
    nPeaks = peakCounts[2];
    xs = newArray(nPeaks);
    ys = newArray(nPeaks);
    vals = newArray(nPeaks);
    for (i = 0; i < nPeaks; i++) {
        xs[i] = xsWork[i];
        ys[i] = ysWork[i];
        vals[i] = valsWork[i];
    }
    if (indexOf(peakMethod, "Masked regions") >= 0) {
        isRegionMethodUsed = 1;
        detectionType = "regions";
        regAreaExport = newArray(nPeaks);
        regMeanExport = newArray(nPeaks);
        regBXExport = newArray(nPeaks);
        regBYExport = newArray(nPeaks);
        regWExport = newArray(nPeaks);
        regHExport = newArray(nPeaks);
        regFeretExport = newArray(nPeaks);
        regMinFeretExport = newArray(nPeaks);
        regCircExport = newArray(nPeaks);
        regDistMaskExport = newArray(nPeaks);
        regDistMaskCentExport = newArray(nPeaks);
        regDistMidExport = newArray(nPeaks);
        for (i = 0; i < nPeaks; i++) {
            regAreaExport[i] = peakRegAreaWork[i];
            regMeanExport[i] = peakRegMeanWork[i];
            regBXExport[i] = peakRegBXWork[i];
            regBYExport[i] = peakRegBYWork[i];
            regWExport[i] = peakRegWWork[i];
            regHExport[i] = peakRegHWork[i];
            regFeretExport[i] = peakRegFeretWork[i];
            regMinFeretExport[i] = peakRegMinFeretWork[i];
            regCircExport[i] = peakRegCircWork[i];
            regDistMaskExport[i] = peakRegDistMaskWork[i];
            regDistMaskCentExport[i] = peakRegDistMaskCentWork[i];
            regDistMidExport[i] = peakRegDistMidWork[i];
        }
    } else {
        isRegionMethodUsed = 0;
        detectionType = "peaks";
        regAreaExport = newArray(0);
        regMeanExport = newArray(0);
        regBXExport = newArray(0);
        regBYExport = newArray(0);
        regWExport = newArray(0);
        regHExport = newArray(0);
        regFeretExport = newArray(0);
        regMinFeretExport = newArray(0);
        regCircExport = newArray(0);
        regDistMaskExport = newArray(0);
        regDistMaskCentExport = newArray(0);
        regDistMidExport = newArray(0);
    }
}

macro "Neurite FluoroQuant - Fast Compute..." {
    print("Neurite FluoroQuant - Fast Compute " + MACRO_VERSION);
    showStatus("Neurite FluoroQuant - Fast Compute " + MACRO_VERSION);
    idAnalysisCorr = 0;
    idAnalysis2Corr = 0;
    idMaskCore = 0;
    hiddenComputeOn = 0;
    // ----- Load previous settings -----
    def_outParent  = prefGetString("unified.outParent", "");
    def_autoIncrementOutput = prefGetBoolean("unified.autoIncrementOutput", true);
    def_fovUm      = prefGetNumber("unified.fovUm", 0);
    def_roiIdx     = round(prefGetNumber("unified.roiIdx", 1));
    def_ch2Idx     = round(prefGetNumber("unified.ch2Idx", 2));
    def_straightUm = prefGetNumber("unified.straightUm", 3.0);
    def_tissueType = prefGetString("unified.tissueType", "2D (single plane)");
    if (def_tissueType == "3D (z-stack)") def_tissueType = "3D (Z-range projection)";
    def_zProjMethod = prefGetString("unified.zProjMethod", "Max Intensity");
    def_compensateZIntensity = prefGetBoolean("unified.compensateZIntensity", true);
    def_rbUm       = prefGetNumber("unified.rbUm", 3.0);
    def_maskThr    = prefGetNumber("unified.maskThr", 0.2);
    def_maskDilUm  = prefGetNumber("unified.maskDilUm", 0.0);
    def_includeCh2IntSummary = prefGetBoolean("unified.includeCh2IntSummary", true);
    def_includeNP  = prefGetBoolean("unified.includeNP", false);
    legacyIncludeCv = prefGetBoolean("unified.includeCV", true);
    def_includeCvSummary = prefGetBoolean("unified.includeCvSummary", legacyIncludeCv);
    def_wcvUm      = prefGetNumber("unified.wcvUm", 10.0);
    def_doDetrend  = prefGetBoolean("unified.doDetrend", true);
    def_detrendUm  = prefGetNumber("unified.detrendUm", 20.0);
    def_includeMaskWidth = prefGetBoolean("unified.includeMaskWidth", true);
    def_includeCh2Profile = prefGetBoolean("unified.includeCh2Profile", true);
    legacyIncludeCh1 = prefGetBoolean("unified.includeCh1", true);
    def_includeRoiSummary = prefGetBoolean("unified.includeRoiSummary", legacyIncludeCh1);
    def_includeRoiProfile = prefGetBoolean("unified.includeRoiProfile", legacyIncludeCh1);
    def_includeAnalysisParams = prefGetBoolean("unified.includeAnalysisParams", true);
    def_peakMethod = "Masked regions (threshold + particles)";
    def_tolPct     = prefGetNumber("unified.tolPct", 12.0);
    def_minSepUm   = prefGetNumber("unified.minSepUm", 0.7);
    def_markerUm   = prefGetNumber("unified.markerUm", 1.5);
    def_doRefine   = prefGetBoolean("unified.doRefine", true);
    def_regionMinAreaUm = prefGetNumber("unified.regionMinAreaUm", 0.05);
    def_regionMaxAreaUm = prefGetNumber("unified.regionMaxAreaUm", 0);
    def_regionSmoothUm  = prefGetNumber("unified.regionSmoothUm", 0);
    def_regionPromPct   = prefGetNumber("unified.regionPromPct", 15.0);
    def_regionDogUm     = prefGetNumber("unified.regionDogUm", 4.0);
    def_batchMode       = prefGetBoolean("unified.batchMode", false);
    def_lastInputDir    = prefGetString("unified.lastInputDir", def_outParent);
    def_lastInputPath   = prefGetString("unified.lastInputPath", "");
    def_usePeakDetection = prefGetBoolean("unified.doPeakDetection", true);
    def_projZFirst      = round(prefGetNumber("unified.projZFirst", 1));
    def_projZLast       = round(prefGetNumber("unified.projZLast", 1));
    def_appendMasterSummary = prefGetBoolean("unified.appendMasterSummary", false);
    def_masterSummaryPath   = prefGetString("unified.masterSummaryPath", "");
    def_adjacentMode        = prefGetBoolean("unified.adjacentMode", false);
    def_adjSecondAnalysisIdx = round(prefGetNumber("unified.adjSecondAnalysisIdx", 0));
    def_adjMaxDistMidUm     = prefGetNumber("unified.adjacentMaxDistMidUm", 1.5);
    // Second-channel detection prefs (adjacent dual pass B only).
    def_adj2_tolPct         = prefGetNumber("unified.adj2.tolPct", def_tolPct);
    def_adj2_minSepUm       = prefGetNumber("unified.adj2.minSepUm", def_minSepUm);
    def_adj2_regionMinAreaUm = prefGetNumber("unified.adj2.regionMinAreaUm", def_regionMinAreaUm);
    def_adj2_regionMaxAreaUm = prefGetNumber("unified.adj2.regionMaxAreaUm", def_regionMaxAreaUm);
    def_adj2_regionSmoothUm  = prefGetNumber("unified.adj2.regionSmoothUm", def_regionSmoothUm);
    def_adj2_regionPromPct   = prefGetNumber("unified.adj2.regionPromPct", def_regionPromPct);
    def_adj2_regionDogUm     = prefGetNumber("unified.adj2.regionDogUm", def_regionDogUm);
    // The rejection distance is edited in the detection tuning loop, not in a setup dialog, so seed
    // it from preferences here for both the interactive and batch paths.
    adjMaxDistMidUm  = def_adjMaxDistMidUm;
    if (adjMaxDistMidUm < 0) adjMaxDistMidUm = 0;

    // Seed the input file field with the exact last file (so accepting it reopens
    // that file), falling back to the last folder only so Browse starts there.
    if (def_lastInputPath != "") def_inputSeed = def_lastInputPath; else def_inputSeed = def_lastInputDir;

    // ----- Parameter dialogs (split for clarity and to prevent UI overflow) -----
    Dialog.create("Neurite FluoroQuant - Fast Compute - Input / Output selection");
    Dialog.addMessage("Choose the input image and output folder.");
    Dialog.addFile("Open 2-channel image:", def_inputSeed);
    Dialog.addDirectory("Output location:", def_outParent);
    Dialog.addCheckbox("Auto-increment output folder", def_autoIncrementOutput);
    Dialog.addMessage("Each run creates a new numbered folder (..._output_001, ..._output_002, ...) when auto-increment is on.");
    Dialog.addMessage(" ");
    Dialog.addCheckbox("Batch mode", def_batchMode);
    Dialog.addMessage("Batch mode reuses saved settings from previous run; you choose the ROI and tune the mask.");
    Dialog.addMessage(" ");
    Dialog.addCheckbox("Append summary to master CSV", def_appendMasterSummary);
    Dialog.addFile("Master summary CSV:", def_masterSummaryPath);
    Dialog.show();

    inPath    = Dialog.getString();
    outParent = Dialog.getString();
    autoIncrementOutput = Dialog.getCheckbox();
    batchMode = Dialog.getCheckbox();
    appendMasterSummary = Dialog.getCheckbox();
    masterSummaryPath = Dialog.getString();

    // Persist the input/output/master-CSV selections immediately so they stay
    // "sticky" from run to run. Previously these were only saved near the end of
    // all the setup dialogs, so cancelling any later dialog (or aborting the run)
    // discarded them and the next run reverted to an earlier run's paths.
    prefSet("unified.batchMode", batchMode);
    prefSet("unified.outParent", outParent);
    prefSet("unified.autoIncrementOutput", autoIncrementOutput);
    prefSet("unified.appendMasterSummary", appendMasterSummary);
    prefSet("unified.masterSummaryPath", masterSummaryPath);
    if (inPath != "") prefSet("unified.lastInputDir", File.getDirectory(inPath));
    if (File.exists(inPath) && !File.isDirectory(inPath)) prefSet("unified.lastInputPath", inPath);

    // Analysis mode: interactive run asks; batch reuses the last-run preference.
    if (!batchMode) {
        if (def_adjacentMode) modeDefault = "Mask-adjacent";
        else modeDefault = "Within-mask";
        Dialog.create("Neurite FluoroQuant - Fast Compute - Analysis Mode");
        Dialog.addChoice("Analysis mode:", newArray("Within-mask", "Mask-adjacent"), modeDefault);
        Dialog.addMessage("Within-mask: quantify signal inside the neurite mask (e.g. axonal fluorescence).");
        Dialog.addMessage("Mask-adjacent: quantify signal adjacent to the neurite mask (e.g. synaptic contacts).");
        Dialog.show();
        analysisModeChoice = Dialog.getChoice();
        adjacentMode = 0;
        if (analysisModeChoice == "Mask-adjacent") adjacentMode = 1;
        prefSet("unified.adjacentMode", adjacentMode);
        def_adjacentMode = adjacentMode;
    } else {
        adjacentMode = def_adjacentMode;
    }

    if (!batchMode) {
    Dialog.create("Neurite FluoroQuant - Fast Compute - General Settings");
    Dialog.addMessage("Set channel roles, scaling, straightening, and ROI selection method.");
    Dialog.addNumber("ROI channel:", def_roiIdx);
    Dialog.addNumber("Analysis channel:", def_ch2Idx);
    if (adjacentMode) Dialog.addNumber("Second analysis channel (0=none):", def_adjSecondAnalysisIdx);
    Dialog.addNumber("Field of view width (um) [0=use calibration]:", def_fovUm);
    Dialog.addNumber("Straighten width (um):", def_straightUm, 1, 5, "");
    Dialog.addChoice("ROI selection method:", newArray("2D (single plane)", "3D (Z-range projection)", "3D (per-plane trace)"), def_tissueType);
    Dialog.addCheckbox("Compensate intensity between Z segments (3D per-plane trace only; disabled for other methods)", def_compensateZIntensity);
    Dialog.show();

    roiIdx    = Dialog.getNumber();
    ch2Idx    = Dialog.getNumber();
    if (adjacentMode) {
        adjSecondAnalysisIdx = round(Dialog.getNumber());
        if (adjSecondAnalysisIdx < 0) adjSecondAnalysisIdx = 0;
        def_adjSecondAnalysisIdx = adjSecondAnalysisIdx;
    } else {
        adjSecondAnalysisIdx = 0;
    }
    fovUm     = Dialog.getNumber();
    straightenUm = Dialog.getNumber();
    tissueType = Dialog.getChoice();
    compensateZIntensityPref = Dialog.getCheckbox();
    // Dual analysis is adjacent-only; keep def_adjSecondAnalysisIdx for prefs, clear runtime if unused.
    adjDualMode = 0;
    if (adjacentMode && adjSecondAnalysisIdx > 0 && adjSecondAnalysisIdx != ch2Idx) adjDualMode = 1;
    if (!adjDualMode) adjSecondAnalysisIdx = 0;
    is3DProject = (indexOf(tissueType, "Z-range") >= 0);
    is3DPerPlane = (indexOf(tissueType, "per-plane") >= 0);
    is3D       = is3DProject || is3DPerPlane;
    if (is3DPerPlane) compensateZIntensity = compensateZIntensityPref;
    else compensateZIntensity = false;
    maskIdx   = roiIdx;
    zProjMethod = def_zProjMethod;
    rbUm         = def_rbUm;
    maskThr      = def_maskThr;
    maskDilUm    = def_maskDilUm;

    if (adjacentMode) {
    // Adjacent-particle mode produces no in-mask intensity, ROI or CV metrics, so only the
    // two exports that still apply are offered. The unasked preferences are carried through
    // unchanged so a later ordinary run finds them exactly as the user last left them.
    Dialog.create("Neurite FluoroQuant - Fast Compute - Output Metrics Settings");
    Dialog.addMessage("Optionally select additional outputs.");
    addMetricCheckbox("Neurite width", def_includeMaskWidth, "Width at each position along the neurite");
    addMetricCheckbox("Export analysis parameters", def_includeAnalysisParams, "Settings and options recorded for this run");
    Dialog.show();

    includeMaskWidth      = Dialog.getCheckbox();
    includeAnalysisParams = Dialog.getCheckbox();
    includeCh2IntSummary  = def_includeCh2IntSummary;
    includeCh2Profile     = def_includeCh2Profile;
    includeNP             = def_includeNP;
    includeCvSummary      = def_includeCvSummary;
    includeRoiSummary     = def_includeRoiSummary;
    includeRoiProfile     = def_includeRoiProfile;
    } else {
    Dialog.create("Neurite FluoroQuant - Fast Compute - Output Metrics Settings");
    Dialog.addMessage("Choose which metrics are written to the exported CSV files.");
    Dialog.addMessage("Analysis channel", 13, "");
    addMetricCheckbox("Analysis channel summary metrics", def_includeCh2IntSummary, "Mean and integrated analysis signal in mask");
    addMetricCheckbox("Analysis channel full profile", def_includeCh2Profile, "Per-position in-mask intensity along arc");
    addMetricCheckbox("Neuropil", def_includeNP, "Background signal outside mask (summary)");
    addMetricCheckbox("Analysis channel intensity CV", def_includeCvSummary, "Global intensity unevenness (rows in analysis summary)");
    Dialog.addMessage("ROI channel", 13, "");
    addMetricCheckbox("ROI channel summary metrics", def_includeRoiSummary, "Length, intensity CV, neurite width CV, and integrated ROI signal");
    addMetricCheckbox("ROI channel full profile", def_includeRoiProfile, "Per-position ROI intensity along arc");
    addMetricCheckbox("Neurite width", def_includeMaskWidth, "Width at each position along the neurite");
    Dialog.addMessage("Parameters", 13, "");
    addMetricCheckbox("Export analysis parameters", def_includeAnalysisParams, "Settings and options recorded for this run");
    Dialog.show();

    includeCh2IntSummary  = Dialog.getCheckbox();
    includeCh2Profile     = Dialog.getCheckbox();
    includeNP             = Dialog.getCheckbox();
    includeCvSummary      = Dialog.getCheckbox();
    includeRoiSummary     = Dialog.getCheckbox();
    includeRoiProfile     = Dialog.getCheckbox();
    includeMaskWidth      = Dialog.getCheckbox();
    includeAnalysisParams = Dialog.getCheckbox();
    }
    includeCv = includeCvSummary;
    wcvUm     = def_wcvUm;
    doDetrend = def_doDetrend;
    detrendUm = def_detrendUm;
    includeCh1  = includeRoiSummary || includeRoiProfile;
    } else {
        applySettingsFromDefs();
    }

    // ----- Save preferences -----
    prefSet("unified.batchMode", batchMode);
    prefSet("unified.outParent", outParent);
    prefSet("unified.autoIncrementOutput", autoIncrementOutput);
    prefSet("unified.fovUm", fovUm);
    prefSet("unified.roiIdx", roiIdx);
    prefSet("unified.ch2Idx", ch2Idx);
    prefSet("unified.maskIdx", roiIdx);
    prefSet("unified.straightUm", straightenUm);
    prefSet("unified.tissueType", tissueType);
    prefSet("unified.zProjMethod", zProjMethod);
    prefSet("unified.compensateZIntensity", compensateZIntensityPref);
    prefSet("unified.maskThr", maskThr);
    prefSet("unified.maskDilUm", maskDilUm);
    prefSet("unified.includeCh2IntSummary", includeCh2IntSummary);
    prefSet("unified.includeCh2Profile", includeCh2Profile);
    prefSet("unified.includeNP", includeNP);
    prefSet("unified.includeCvSummary", includeCvSummary);
    prefSet("unified.includeCV", includeCv);
    prefSet("unified.includeRoiSummary", includeRoiSummary);
    prefSet("unified.includeRoiProfile", includeRoiProfile);
    prefSet("unified.includeMaskWidth", includeMaskWidth);
    prefSet("unified.includeAnalysisParams", includeAnalysisParams);
    prefSet("unified.includeCh1", includeCh1);
    prefSet("unified.appendMasterSummary", appendMasterSummary);
    prefSet("unified.masterSummaryPath", masterSummaryPath);
    prefSet("unified.adjacentMode", adjacentMode);
    prefSet("unified.adjSecondAnalysisIdx", def_adjSecondAnalysisIdx);

    // Adjacent-particle mode measures nothing inside the mask, so every in-mask metric is
    // switched off here — after the preferences above have been written, so the user's own
    // choices survive untouched for the next ordinary run.
    if (adjacentMode) {
        includeCh2IntSummary = 0;
        includeCh2Profile = 0;
        includeNP = 0;
        includeCvSummary = 0;
        includeCv = 0;
        includeRoiSummary = 0;
        includeRoiProfile = 0;
        includeCh1 = 0;
    }

    batchProcessedCount = 0;
    batchMore = true;
    while (batchMore) {
    if (inPath == "") exit("No image selected.");
    // Guard against a folder or missing path reaching open(): opening a directory
    // imports every image in it as one merged Z-stack. Reject it with a clear message.
    if (File.isDirectory(inPath)) {
        exit("The selected input is a folder, not an image file:\n" + inPath +
             "\nPlease choose the .lif or image file (use Browse).");
    }
    if (!File.exists(inPath)) {
        exit("The selected input does not exist:\n" + inPath +
             "\nPlease choose a valid .lif or image file (use Browse).");
    }
    segZ = newArray(0);
    pathSegCount = 0;
    pathTotalPx = 0;

    // ----- Open and split channels -----
    open(inPath);
    prefSet("unified.lastInputDir", File.getDirectory(inPath));
    if (File.exists(inPath) && !File.isDirectory(inPath)) prefSet("unified.lastInputPath", inPath);
    origTitle = getTitle();

    // ----- Prepare output directory -----
    // Derived after opening so multi-image containers (e.g. Leica .lif) can adopt
    // the individual stack/series name from the window title; single-image files
    // fall back to the input file name exactly as before.
    base = deriveOutputBaseName(inPath, origTitle);
    if (outParent == "") outParent = File.getDirectory(inPath);
    if (!endsWith(outParent, File.separator)) outParent = outParent + File.separator;
    if (autoIncrementOutput) {
        outDir = allocateNumberedOutputDir(outParent, base);
    } else {
        outDir = outParent + base + "_output";
        ensureDirectoryExists(outDir);
    }
    activeOutDir = outDir;
    getDimensions(w, h, ch, sl, fr);
    getPixelSize(unit, calW, calH, calD);
    useCal = (fovUm == 0 && calW > 0);
    if (useCal) umPerPx = umPerPixelFromCalibration(); else umPerPx = fovUm / w;

    // Validate tissue type against the actual Z depth.
    if (is3D && sl <= 1) {
        cancelRunAndExit(outDir, "3D tissue selected, but this image has no Z slices (single plane).\n" +
             "Open a z-stack/hyperstack, or choose 2D tissue.");
    }
    if (!is3D && sl > 1) {
        showMessage("2D tissue on a z-stack",
            "2D tissue selected, but this image has " + sl + " Z slices.\n" +
            "Only the current/first plane will be used.");
    }

    if (ch > 1) {
        run("Split Channels");
    } else {
        // Single-channel image: nothing to split. Both ROI and analysis
        // functions necessarily operate on the sole channel, so map every
        // role to channel 1 and give the window the C1- prefix the rest of
        // the pipeline expects.
        roiIdx = 1;
        ch2Idx = 1;
        maskIdx = 1;
        adjSecondAnalysisIdx = 0;
        adjDualMode = 0;
        rename("C1-" + origTitle);
    }
    /*
     * Determine the titles of the split channel images in a more
     * forgiving manner.  ImageJ names split channels starting
     * with "C<index>-" followed by the base title.  However the
     * exact base title may differ (extension removed, truncated,
     * duplicate counters, etc.), so relying on the original title
     * string can cause failures when selecting the window.  We
     * simply pick the first open window whose title starts with
     * the appropriate channel index.  This makes the macro robust
     * against variations in how ImageJ constructs the channel
     * names.
     */
    roiTitleFound = findSplitChannelTitle(roiIdx);
    ch2TitleFound = findSplitChannelTitle(ch2Idx);
    maskTitleFound = findSplitChannelTitle(maskIdx);
    adj2TitleFound = "";
    if (adjDualMode) adj2TitleFound = findSplitChannelTitle(adjSecondAnalysisIdx);
    /*
     * If for some reason no window starting with the desired channel
     * number exists (e.g. the image is not multi‑channel), fall back
     * to the default naming convention.  We do not require the
     * original title string here, to avoid mismatches with file
     * extensions or truncated names.
     */
    if (roiTitleFound == "") roiTitleFound = "C" + d2s(roiIdx, 0) + "-";
    if (ch2TitleFound == "") ch2TitleFound = "C" + d2s(ch2Idx, 0) + "-";
    if (maskTitleFound == "") maskTitleFound = "C" + d2s(maskIdx, 0) + "-";
    if (adjDualMode && adj2TitleFound == "") adj2TitleFound = "C" + d2s(adjSecondAnalysisIdx, 0) + "-";
    if (!isOpen(roiTitleFound)) {
        cancelRunAndExit(outDir, "Could not find ROI channel window starting with C" + roiIdx + "\n" +
             "Please verify that the image has at least " + roiIdx + " channels and that the channel assignment is correct.");
    }
    if (!isOpen(ch2TitleFound)) {
        cancelRunAndExit(outDir, "Could not find analysis channel window starting with C" + ch2Idx + "\n" +
             "Please verify that the image has at least " + ch2Idx + " channels and that the channel assignment is correct.");
    }
    if (!isOpen(maskTitleFound)) {
        cancelRunAndExit(outDir, "Could not find mask channel window starting with C" + maskIdx + "\n" +
             "Please verify that the image has at least " + maskIdx + " channels and that the channel assignment is correct.");
    }
    if (adjDualMode && !isOpen(adj2TitleFound)) {
        cancelRunAndExit(outDir, "Could not find second analysis channel window starting with C" + adjSecondAnalysisIdx + "\n" +
             "Please verify that the image has at least " + adjSecondAnalysisIdx + " channels and that the channel assignment is correct.");
    }
    // Rename split channels when possible; duplicate only if the same channel index serves multiple roles.
    maxChIdx = roiIdx;
    if (ch2Idx > maxChIdx) maxChIdx = ch2Idx;
    if (maskIdx > maxChIdx) maxChIdx = maskIdx;
    if (adjDualMode && adjSecondAnalysisIdx > maxChIdx) maxChIdx = adjSecondAnalysisIdx;
    consumedTitles = newArray(maxChIdx + 1);
    for (i = 0; i < consumedTitles.length; i++) consumedTitles[i] = "";
    roiWin = "ROI_C" + d2s(roiIdx, 0);
    analysisWin = "Analysis_C" + d2s(ch2Idx, 0);
    maskChannelWin = "Mask_C" + d2s(maskIdx, 0);
    analysis2Win = "";
    analysis2StraightWin = "";
    analysis2CorrWin = "";
    analysis2PlainWin = "";
    analysis2OverlayWin = "";
    analysis2MaskPreviewWin = "";
    analysis2BgsubPreviewWin = "";
    analysis2ExportTag = "";
    if (adjDualMode) {
        analysis2Win = "Analysis2_C" + d2s(adjSecondAnalysisIdx, 0);
        analysis2StraightWin = analysis2Win + "_straight";
        analysis2CorrWin = analysis2Win + "_corr";
        analysis2PlainWin = analysis2Win + "_plain";
        analysis2OverlayWin = analysis2Win + "_overlay";
        analysis2MaskPreviewWin = analysis2Win + "_mask_preview";
        analysis2BgsubPreviewWin = analysis2Win + "_bgsub_preview";
        analysis2ExportTag = "_analysis_C" + d2s(adjSecondAnalysisIdx, 0);
    }
    roiStraightWin = roiWin + "_straight";
    analysisStraightWin = analysisWin + "_straight";
    roiCorrWin = roiWin + "_corr";
    analysisCorrWin = analysisWin + "_corr";
    analysisNormWin = analysisWin + "_norm";
    roiPlainWin = roiWin + "_plain";
    analysisPlainWin = analysisWin + "_plain";
    analysisOverlayWin = analysisWin + "_overlay";
    maskPreviewWin = maskChannelWin + "_preview";
    analysisMaskPreviewWin = analysisWin + "_mask_preview";
    analysisBgsubPreviewWin = analysisWin + "_bgsub_preview";
    roiBgsubPreviewWin = roiWin + "_bgsub_preview";
    peakPreviewPrefix = analysisWin + "_peak_preview";
    peakProfilePlotWin = analysisWin + "_peak_profile";
    assignRoleChannel(roiIdx, roiTitleFound, roiWin, consumedTitles);
    assignRoleChannel(ch2Idx, ch2TitleFound, analysisWin, consumedTitles);
    assignRoleChannel(maskIdx, maskTitleFound, maskChannelWin, consumedTitles);
    if (adjDualMode) assignRoleChannel(adjSecondAnalysisIdx, adj2TitleFound, analysis2Win, consumedTitles);
    closeLeftoverChannelWindows(origTitle, ch);
    roiExportTag = "_roi_C" + d2s(roiIdx, 0);
    analysisExportTag = "_analysis_C" + d2s(ch2Idx, 0);
    maskExportTag = "_mask_C" + d2s(maskIdx, 0);
    analysisOutDir = ensureOutputSubdir(outDir, "Analysis");
    paramsOutDir = ensureOutputSubdir(outDir, "Parameters");
    roiOutDir = ensureOutputSubdir(outDir, "ROI");
    analysisProfileRelPath = "Analysis" + File.separator + base + analysisExportTag + "_profile.csv";
    if (!batchMode) {
        chAssignMsg = "ROI channel = C" + d2s(roiIdx, 0) + "\n" +
            "Analysis channel = C" + d2s(ch2Idx, 0);
        if (adjDualMode) chAssignMsg = chAssignMsg + "\nSecond analysis channel = C" + d2s(adjSecondAnalysisIdx, 0);
        chAssignMsg = chAssignMsg + "\nMask generated from ROI channel";
        showMessage("Channel assignment", chAssignMsg);
    }

    // ----- ROI tracing on ROI channel -----
    selectWindow(roiWin);
    run("Enhance Contrast", "saturated=0.35");
    straightPx = round(straightenUm / umPerPx); if (straightPx < 2) straightPx = 2;
    if (is3DProject) {
        pickZRangeGlobal(roiWin);
        prefSet("unified.zProjMethod", zProjMethod);
        projectChannelStack(roiWin, projZFirst, projZLast, zProjMethod);
        projectChannelStack(analysisWin, projZFirst, projZLast, zProjMethod);
        projectChannelStack(maskChannelWin, projZFirst, projZLast, zProjMethod);
        if (adjDualMode) projectChannelStack(analysis2Win, projZFirst, projZLast, zProjMethod);
    }
    if (is3DPerPlane) {
        // 3D per-plane: trace as ordered segments (ROI Manager + segZ[]), then align and stitch.
        nSeg3D = capturePathSegments(roiWin);
        roiLenUm = pathTotalPx * umPerPx;
        straightenAlongPath(analysisWin, analysisStraightWin, straightPx, nSeg3D);
        closeWindowIfOpen(analysisWin);
        if (adjDualMode) {
            straightenAlongPath(analysis2Win, analysis2StraightWin, straightPx, nSeg3D);
            closeWindowIfOpen(analysis2Win);
        }
        straightenAlongPath(maskChannelWin, "Mask_straight", straightPx, nSeg3D);
        closeWindowIfOpen(maskChannelWin);
    } else {
        waitForUser("Draw a Segmented Line / Polyline along the structure to analyse, then click OK.");
        rt = selectionType();
        if (rt != 5 && rt != 6) cancelRunAndExit(outDir, "");
        getSelectionCoordinates(px, py);
        if (px.length < 2) cancelRunAndExit(outDir, "");
        total_px = 0;
        for (i = 1; i < px.length; i++) {
            dx = px[i] - px[i-1]; dy = py[i] - py[i-1]; total_px += sqrt(dx*dx + dy*dy);
        }
        roiLenUm = total_px * umPerPx;

        // ----- Straighten analysis channel and mask (ROI deferred until after mask tuning) -----
        roiManager("Reset"); roiManager("Add");
        selectWindow(analysisWin); roiManager("Select", 0); run("Straighten...", "line=" + straightPx); rename(analysisStraightWin);
        closeWindowIfOpen(analysisWin);
        if (adjDualMode) {
            selectWindow(analysis2Win); roiManager("Select", 0); run("Straighten...", "line=" + straightPx); rename(analysis2StraightWin);
            closeWindowIfOpen(analysis2Win);
        }
        selectWindow(maskChannelWin); roiManager("Select", 0); run("Straighten...", "line=" + straightPx); rename("Mask_straight");
        closeWindowIfOpen(maskChannelWin);
    }

    // ----- Background subtraction tuning (before mask; same rbUm used for mask build and metrics) -----
    if (batchMode) {
        rbUm = def_rbUm;
        if (rbUm < 0) rbUm = 0;
    } else {
    doneBgsub = false;
    bgsubPreviewReady = 0;
    setBatchMode(false);
    while (!doneBgsub) {
        Dialog.createNonBlocking("Background subtraction tuning");
        if (adjDualMode) {
            if (bgsubPreviewReady == 0) {
                Dialog.addMessage("Tune rolling-ball background subtraction radius.\n\n" +
                    "Left = ROI (mask source). Middle = primary analysis. Right = second analysis.");
            } else {
                Dialog.addMessage("Left: ROI after bgsub. Middle: primary analysis after bgsub. Right: second analysis after bgsub.");
            }
        } else if (bgsubPreviewReady == 0) {
            Dialog.addMessage("Tune rolling-ball background subtraction radius.\n\n" +
                "Left = ROI (mask source). Right = analysis (intensity metrics).");
        } else {
            Dialog.addMessage("Left: ROI channel after background subtraction. Right: analysis channel after background subtraction.");
        }
        Dialog.addNumber("Rolling-ball radius (um):", rbUm, 1, 5, "0 = no background subtraction");
        Dialog.addChoice("Action:", newArray("Update preview", "Accept and continue", "Cancel run (discard image)"), "Update preview");
        Dialog.show();
        rbUm = Dialog.getNumber();
        if (rbUm < 0) rbUm = 0;
        actionBgsub = Dialog.getChoice();
        if (indexOf(actionBgsub, "Cancel run") >= 0) cancelRunAndExit(outDir, "");
        doUpdateBgsub = indexOf(actionBgsub, "Update") >= 0;
        doAcceptBgsub = indexOf(actionBgsub, "Accept") >= 0;
        if (doUpdateBgsub) {
            buildBgsubPreview("Mask_straight", roiBgsubPreviewWin, rbUm);
            buildBgsubPreview(analysisStraightWin, analysisBgsubPreviewWin, rbUm);
            if (adjDualMode && analysis2StraightWin != "" && isOpen(analysis2StraightWin)) {
                buildBgsubPreview(analysis2StraightWin, analysis2BgsubPreviewWin, rbUm);
            }
            layoutBgsubPreviews();
            bgsubPreviewReady = 1;
        }
        if (doAcceptBgsub) {
            closeWindowIfOpen(roiBgsubPreviewWin);
            closeWindowIfOpen(analysisBgsubPreviewWin);
            if (analysis2BgsubPreviewWin != "") closeWindowIfOpen(analysis2BgsubPreviewWin);
            prefSet("unified.rbUm", rbUm);
            doneBgsub = true;
        }
    }
    }

    // ----- Generate binary mask on Mask_straight with tuning loop -----
    // Batch mode still tunes the mask per image (threshold varies with quality);
    // rolling-ball radius stays at the saved value from the previous run.
    maskThrFrac = maskThr;
    if (maskThrFrac < 0) maskThrFrac = 0;
    if (maskThrFrac > 1) maskThrFrac = 1;
    doneMask = false;
    while (!doneMask) {
        closeWindowIfOpen(maskPreviewWin);
        closeWindowIfOpen(analysisMaskPreviewWin);
        if (analysis2MaskPreviewWin != "") closeWindowIfOpen(analysis2MaskPreviewWin);

        buildMaskWorkFromStraight(maskThrFrac);
        buildMaskPreview("Mask_straight", "Mask_work", maskPreviewWin);
        buildMaskPreview(analysisStraightWin, "Mask_work", analysisMaskPreviewWin);
        if (adjDualMode && analysis2StraightWin != "" && isOpen(analysis2StraightWin)) {
            buildMaskPreview(analysis2StraightWin, "Mask_work", analysis2MaskPreviewWin);
        }
        closeWindowIfOpen("Mask_work");
        layoutTuningPreviews();

        Dialog.createNonBlocking("Mask threshold tuning");
        Dialog.addMessage("Adjust the threshold until the red mask outline fits the neurite, then accept.");
        Dialog.addNumber("Mask threshold fraction (0-1):", maskThrFrac, 3, 3, "");
        Dialog.addChoice("Action:", newArray("Update preview", "Accept mask and continue", "Cancel run (discard image)"), "Update preview");
        // The built-in "Cancel" button aborts the macro uncatchably (ImageJ default),
        // so a catchable "Cancel run" action is provided to close all windows cleanly.
        Dialog.show();
        newThr = Dialog.getNumber();
        action = Dialog.getChoice();
        if (indexOf(action, "Cancel run") >= 0) cancelRunAndExit(outDir, "");
        if (newThr < 0) newThr = 0;
        if (newThr > 1) newThr = 1;
        maskThrFrac = newThr;
        prefSet("unified.maskThr", maskThrFrac);

        acceptMask = indexOf(action, "Accept") >= 0;
        if (acceptMask) {
            closeWindowIfOpen(maskPreviewWin);
            closeWindowIfOpen(analysisMaskPreviewWin);
            if (analysis2MaskPreviewWin != "") closeWindowIfOpen(analysis2MaskPreviewWin);
            doneMask = true;
        } else {
            closeWindowIfOpen(maskPreviewWin);
            closeWindowIfOpen(analysisMaskPreviewWin);
            if (analysis2MaskPreviewWin != "") closeWindowIfOpen(analysis2MaskPreviewWin);
        }
    }
    maskThr = maskThrFrac;
    if (!isOpen("Mask_work")) buildMaskWorkFromStraight(maskThrFrac);
    // Finalize the mask: rename the working mask to Mask_final
    maskTitle = "Mask_final";
    if (isOpen(maskTitle)) { selectWindow(maskTitle); run("Close"); }
    if (isOpen("Mask_work")) {
        selectWindow("Mask_work");
        rename(maskTitle);
    } else {
        cancelRunAndExit(activeOutDir, "Unable to find the working mask window to rename to Mask_final.\n" +
             "Please verify the mask thresholding step.");
    }
    requireWindow(maskTitle);
    setOption("Changes", false);

    requireWindow(maskTitle);
    beginHiddenCompute();

    // ----- Allocate mask geometry arrays (filled after Mask_core is built) -----
    getDimensions(swM, shM, cM, sM, fM);
    widthsPx = newArray(swM);
    widthsUm = newArray(swM);
    maskFirstYPx = newArray(swM);
    maskLastYPx = newArray(swM);
    maskCorePresent = newArray(swM);
    maskMidYPx = (shM - 1) / 2.0;
    for (x = 0; x < swM; x++) {
        widthsPx[x] = 0; widthsUm[x] = 0;
        maskFirstYPx[x] = -1; maskLastYPx[x] = -1;
        maskCorePresent[x] = 0;
    }

    // Straighten ROI channel when ROI metrics are requested, or when dual-adjacent needs it
    // for the merged RGB overlay (ROI quantitative exports stay off in adjacent mode).
    if (includeCh1 || adjDualMode) {
        if (is3DPerPlane) {
            // Reuse the captured per-plane segments (still in the ROI Manager) so the ROI channel
            // aligns column-for-column with the analysis/mask straightened images.
            straightenAlongPath(roiWin, roiStraightWin, straightPx, pathSegCount);
            closeWindowIfOpen(roiWin);
        } else {
            selectWindow(roiWin);
            if (rt == 5) makeSelection("polyline", px, py);
            else makeSelection("freeline", px, py);
            roiManager("Reset");
            roiManager("Add");
            roiManager("Select", 0);
            run("Straighten...", "line=" + straightPx);
            rename(roiStraightWin);
            closeWindowIfOpen(roiWin);
        }
    } else {
        closeWindowIfOpen(roiWin);
    }

    // ----- Background subtraction on analysis channel (always) and ROI (optional) -----
    activateImage(analysisStraightWin);
    duplicateCurrentAs(analysisCorrWin);
    idAnalysisCorr = pinCurrentImage();
    if (rbUm > 0) {
        rbPx2 = round(rbUm / umPerPx); if (rbPx2 < 1) rbPx2 = 1;
        run("Subtract Background...", "rolling=" + rbPx2 + " sliding");
    }
    if (adjDualMode) {
        activateImage(analysis2StraightWin);
        duplicateCurrentAs(analysis2CorrWin);
        idAnalysis2Corr = pinCurrentImage();
        if (rbUm > 0) {
            rbPx2b = round(rbUm / umPerPx); if (rbPx2b < 1) rbPx2b = 1;
            run("Subtract Background...", "rolling=" + rbPx2b + " sliding");
        }
    }

    if (includeCh1) {
        activateImage(roiStraightWin);
        duplicateCurrentAs(roiCorrWin);
        pinCurrentImage();
        if (rbUm > 0) {
            rbPx1 = round(rbUm / umPerPx); if (rbPx1 < 1) rbPx1 = 1;
            run("Subtract Background...", "rolling=" + rbPx1 + " sliding");
        }
    }
    // Keep a dedicated ROI grayscale for the dual merged overlay after pre-peak cleanup.
    adjMergeRoiWin = "";
    adjRoiZipPass1 = "";
    adjRoiZipPass2 = "";
    if (adjDualMode && imageIsOpen(roiStraightWin)) {
        adjMergeRoiWin = "AdjMerge_ROI";
        closeWindowIfOpen(adjMergeRoiWin);
        if (includeCh1 && imageIsOpen(roiCorrWin)) {
            activateImage(roiCorrWin);
            duplicateCurrentAs(adjMergeRoiWin);
        } else {
            activateImage(roiStraightWin);
            duplicateCurrentAs(adjMergeRoiWin);
            if (rbUm > 0) {
                rbPxMrg = round(rbUm / umPerPx); if (rbPxMrg < 1) rbPxMrg = 1;
                run("Subtract Background...", "rolling=" + rbPxMrg + " sliding");
            }
        }
    }

    // Reduce the mask to the midline neurite; Mask_final is then redundant.
    adjCoreMaskTitle = maskTitle;
    adjUmPerPx = umPerPx;
    buildCoreMask();
    if (adjCoreMaskTitle == "Mask_core" && imageIsOpen("Mask_core")) {
        closeWindowIfOpen("Mask_final");
        maskTitle = "Mask_core";
        adjCoreMaskTitle = "Mask_core";
    }
    // Sole mask export: the core neurite used for all downstream metrics.
    requireWindow(maskTitle);
    saveAs("PNG", roiOutDir + File.separator + base + maskExportTag + "_mask.png");
    if (indexOf(getTitle(), maskTitle) != 0) rename(maskTitle);

    // Widths from the core mask only (maskCorePresent excludes gap-filled empty columns).
    sumW = 0; cntW = 0;
    for (x = 0; x < swM; x++) {
        if (maskCorePresent[x] == 1 && maskFirstYPx[x] >= 0) {
            wpx = maskLastYPx[x] - maskFirstYPx[x] + 1;
            widthsPx[x] = wpx;
            wum = wpx * umPerPx;
            widthsUm[x] = wum;
            sumW += wum;
            cntW++;
        } else {
            widthsPx[x] = 0; widthsUm[x] = 0;
        }
    }
    if (cntW > 0) maskWidthMeanUm = sumW / cntW; else maskWidthMeanUm = 0;
    if (includeRoiSummary) maskWidthCv = cvOfValidValues(widthsUm, maskCorePresent, widthsUm.length);
    else maskWidthCv = 0;
    roiCvGlobal = 0;

    // Save width table (optional export)
    if (includeMaskWidth) {
        Table.create("Neurite Width");
        for (x = 0; x < swM; x++) {
            if (swM > 1) xnorm = x / (swM - 1); else xnorm = 0;
            xum = xnorm * roiLenUm;
            Table.set("X px straight", x, x);
            Table.set("X arc (um)", x, xum);
            Table.set("X norm 0-1", x, xnorm);
            Table.set("Width px", x, widthsPx[x]);
            Table.set("Width (um)", x, widthsUm[x]);
        }
        Table.save(roiOutDir + File.separator + base + "_widths.csv");
        if (isOpen("Neurite Width")) { selectWindow("Neurite Width"); run("Close"); }
    }

    // Adjacent-particle mode quantifies signal outside the mask, so every in-mask intensity,
    // ROI and CV measurement below is skipped; the run continues straight to patch detection.
    if (!adjacentMode) {

    // ----- Measure Ch2 intensity metrics inside mask -----
    // Compute raw mean inside mask and neuropil per X
    activateImage(analysisStraightWin); getDimensions(swC2, shC2, cC2, sC2, fC2);
    requireWindow(maskTitle); getDimensions(swM2, shM2, cM2, sM2, fM2);
    // Ensure same width
    swUse = swC2; if (swM2 < swUse) swUse = swM2;
    shUse = shC2; if (shM2 < shUse) shUse = shM2;
    ch2_inMask_raw  = newArray(swUse);
    ch2_neuropil_raw = newArray(swUse);
    ch2_corr_raw    = newArray(swUse);
    ch2_enrich_raw  = newArray(swUse);
    dummyValidRaw = newArray(swUse);
    fillMaskedColumnProfiles(analysisStraightWin, maskTitle, swUse, shUse, ch2_inMask_raw, ch2_neuropil_raw, dummyValidRaw);
    ch2_sum_inMask_raw = scanSumIn; ch2_cnt_inMask_raw = scanCntIn;
    ch2_sum_neuropil_raw = scanSumOut; ch2_cnt_neuropil_raw = scanCntOut;
    if (includeNP) {
        for (x = 0; x < swUse; x++) {
            ch2_corr_raw[x] = ch2_inMask_raw[x] - ch2_neuropil_raw[x];
            if (ch2_neuropil_raw[x] > 0) ch2_enrich_raw[x] = ch2_inMask_raw[x] / ch2_neuropil_raw[x]; else ch2_enrich_raw[x] = 0;
        }
    }
    if (ch2_cnt_inMask_raw > 0) ch2_mean_inMask_raw = ch2_sum_inMask_raw / ch2_cnt_inMask_raw; else ch2_mean_inMask_raw = 0;
    if (ch2_cnt_neuropil_raw > 0) ch2_mean_neuropil_raw = ch2_sum_neuropil_raw / ch2_cnt_neuropil_raw; else ch2_mean_neuropil_raw = 0;
    if (includeNP) {
        ch2_mean_corrected_raw = ch2_mean_inMask_raw - ch2_mean_neuropil_raw;
        if (ch2_mean_neuropil_raw > 0) ch2_mean_enrich_raw = ch2_mean_inMask_raw / ch2_mean_neuropil_raw; else ch2_mean_enrich_raw = 0;
    }

    // ----- Measure Ch2 intensity metrics on background-subtracted image inside mask -----
    activateImage(analysisCorrWin); getDimensions(swC2b, shC2b, cC2b, sC2b, fC2b);
    swUse2 = swC2b; if (swM2 < swUse2) swUse2 = swM2;
    shUse2 = shC2b; if (shM2 < shUse2) shUse2 = shM2;
    ch2_inMask_bgsub  = newArray(swUse2);
    ch2_inMask_valid  = newArray(swUse2);
    ch2_neuropil_bgsub = newArray(swUse2);
    ch2_corr_bgsub    = newArray(swUse2);
    ch2_enrich_bgsub  = newArray(swUse2);
    fillMaskedColumnProfiles(analysisCorrWin, maskTitle, swUse2, shUse2, ch2_inMask_bgsub, ch2_neuropil_bgsub, ch2_inMask_valid);
    ch2_sum_inMask_bgsub = scanSumIn; ch2_cnt_inMask_bgsub = scanCntIn;
    ch2_sum_neuropil_bgsub = scanSumOut; ch2_cnt_neuropil_bgsub = scanCntOut;
    if (includeNP) {
        for (x = 0; x < swUse2; x++) {
            ch2_corr_bgsub[x] = ch2_inMask_bgsub[x] - ch2_neuropil_bgsub[x];
            if (ch2_neuropil_bgsub[x] > 0) ch2_enrich_bgsub[x] = ch2_inMask_bgsub[x] / ch2_neuropil_bgsub[x]; else ch2_enrich_bgsub[x] = 0;
        }
    }
    if (ch2_cnt_inMask_bgsub > 0) ch2_mean_inMask_bgsub = ch2_sum_inMask_bgsub / ch2_cnt_inMask_bgsub; else ch2_mean_inMask_bgsub = 0;
    if (ch2_cnt_neuropil_bgsub > 0) ch2_mean_neuropil_bgsub = ch2_sum_neuropil_bgsub / ch2_cnt_neuropil_bgsub; else ch2_mean_neuropil_bgsub = 0;
    if (includeNP) {
        ch2_mean_corrected_bgsub = ch2_mean_inMask_bgsub - ch2_mean_neuropil_bgsub;
        if (ch2_mean_neuropil_bgsub > 0) ch2_mean_enrich_bgsub = ch2_mean_inMask_bgsub / ch2_mean_neuropil_bgsub; else ch2_mean_enrich_bgsub = 0;
    }

    // ----- Measure ROI channel intensity metrics inside mask -----
    if (includeCh1) {
        if (roiIdx == ch2Idx) {
            profileC1Raw = newArray(swUse);
            for (i = 0; i < swUse; i++) profileC1Raw[i] = ch2_inMask_raw[i];
            profileC1Bg = newArray(swUse2);
            for (i = 0; i < swUse2; i++) profileC1Bg[i] = ch2_inMask_bgsub[i];
            meanC1Raw = ch2_mean_inMask_raw;
            meanC1Bg = ch2_mean_inMask_bgsub;
            projIntC1Raw = ch2_sum_inMask_raw;
            projIntC1Bg = ch2_sum_inMask_bgsub;
            activateImage(roiStraightWin); getDimensions(swC1, shC1, cC1, sC1, fC1);
            nptsC1 = swUse;
        } else {
            activateImage(roiStraightWin); getDimensions(swC1, shC1, cC1, sC1, fC1);
            swUseC1 = swC1; if (swM2 < swUseC1) swUseC1 = swM2;
            shUseC1 = shC1; if (shM2 < shUseC1) shUseC1 = shM2;
            profileC1Raw = newArray(swUseC1);
            dummyNpC1Raw = newArray(swUseC1);
            dummyValidC1Raw = newArray(swUseC1);
            fillMaskedColumnProfiles(roiStraightWin, maskTitle, swUseC1, shUseC1, profileC1Raw, dummyNpC1Raw, dummyValidC1Raw);
            c1_sum_inMask_raw = scanSumIn; c1_cnt_inMask_raw = scanCntIn;
            if (c1_cnt_inMask_raw > 0) meanC1Raw = c1_sum_inMask_raw / c1_cnt_inMask_raw; else meanC1Raw = 0;
            projIntC1Raw = c1_sum_inMask_raw;
            activateImage(roiCorrWin); getDimensions(swC1c, shC1c, cC1c, sC1c, fC1c);
            swUseC1b = swC1c; if (swM2 < swUseC1b) swUseC1b = swM2;
            shUseC1b = shC1c; if (shM2 < shUseC1b) shUseC1b = shM2;
            profileC1Bg = newArray(swUseC1b);
            c1_inMask_valid = newArray(swUseC1b);
            dummyNpC1Bg = newArray(swUseC1b);
            fillMaskedColumnProfiles(roiCorrWin, maskTitle, swUseC1b, shUseC1b, profileC1Bg, dummyNpC1Bg, c1_inMask_valid);
            c1_sum_inMask_bgsub = scanSumIn; c1_cnt_inMask_bgsub = scanCntIn;
            if (c1_cnt_inMask_bgsub > 0) meanC1Bg = c1_sum_inMask_bgsub / c1_cnt_inMask_bgsub; else meanC1Bg = 0;
            projIntC1Bg = c1_sum_inMask_bgsub;
            nptsC1 = swUseC1;
        }
        if (includeRoiSummary) {
            if (roiIdx == ch2Idx) roiCvGlobal = cvOfValidValues(profileC1Bg, ch2_inMask_valid, profileC1Bg.length);
            else roiCvGlobal = cvOfValidValues(profileC1Bg, c1_inMask_valid, profileC1Bg.length);
        } else {
            roiCvGlobal = 0;
        }
    } else {
        roiCvGlobal = 0;
    }

    // Mask-restricted projection integrated intensity (sum of in-mask pixel values)
    projIntC2Bg  = ch2_sum_inMask_bgsub;
    projIntC2Raw = ch2_sum_inMask_raw;
    swC2c = swUse2;

    // ----- Global CV on in-mask bgsub column means -----
    // Uses the same per-column series as Analysis Channel Intensity Bgsub. Empty mask
    // columns are omitted from global / detrended statistics.
    if (includeCv) {
        nptsN = ch2_inMask_bgsub.length;
        sumN = 0; sumsqN = 0; nValidCv = 0;
        for (i = 0; i < nptsN; i++) {
            if (ch2_inMask_valid[i]) {
                v = ch2_inMask_bgsub[i];
                sumN += v;
                sumsqN += v * v;
                nValidCv++;
            }
        }
        if (nValidCv > 0) meanCv = sumN / nValidCv; else meanCv = 0;
        if (nValidCv > 1) {
            varCv = (sumsqN / nValidCv) - (meanCv * meanCv);
            if (varCv < 0) varCv = 0;
            sdCv = sqrt(varCv);
        } else {
            sdCv = 0;
        }
        if (meanCv > 0) cv_global = sdCv / meanCv; else cv_global = 0;
        if (doDetrend) {
            dw = round(detrendUm / umPerPx); if (dw < 3) dw = 3; if ((dw % 2) == 0) dw++;
            dh = floor(dw / 2);
            runmean = newArray(nptsN);
            prof_detr = newArray(nptsN);
            for (i = 0; i < nptsN; i++) {
                runmean[i] = 0;
                prof_detr[i] = 0;
                if (ch2_inMask_valid[i]) {
                    s = i - dh; if (s < 0) s = 0;
                    e = i + dh; if (e >= nptsN) e = nptsN - 1;
                    cnt = 0; msum = 0;
                    for (k = s; k <= e; k++) {
                        if (ch2_inMask_valid[k]) { msum += ch2_inMask_bgsub[k]; cnt++; }
                    }
                    if (cnt > 0) {
                        m = msum / cnt;
                        runmean[i] = m;
                        prof_detr[i] = ch2_inMask_bgsub[i] - m;
                    }
                }
            }
            sumd = 0; sumd2 = 0; nValidD = 0;
            for (i = 0; i < nptsN; i++) {
                if (ch2_inMask_valid[i]) {
                    sumd += prof_detr[i];
                    sumd2 += prof_detr[i] * prof_detr[i];
                    nValidD++;
                }
            }
            if (nValidD > 1) {
                meand = sumd / nValidD;
                vard = (sumd2 / nValidD) - (meand * meand);
                if (vard < 0) vard = 0;
                stdd = sqrt(vard);
            } else {
                stdd = 0;
            }
            if (meanCv > 0) cv_detrended = stdd / meanCv; else cv_detrended = 0;
        }
    }
    closeWindowIfOpen(analysisNormWin);

    }

    saveStraightenedPlainPngs();
    closePrePeakTuningWindows();
    endHiddenCompute();

    // ----- Optional interactive patch detection (after all other measurements) -----
    regionExportRoiPath = "";
    xs = newArray(0);
    ys = newArray(0);
    vals = newArray(0);
    regAreaExport = newArray(0);
    regMeanExport = newArray(0);
    regBXExport = newArray(0);
    regBYExport = newArray(0);
    regWExport = newArray(0);
    regHExport = newArray(0);
    regFeretExport = newArray(0);
    regMinFeretExport = newArray(0);
    regCircExport = newArray(0);
    regDistMaskExport = newArray(0);
    regDistMaskCentExport = newArray(0);
    regDistMidExport = newArray(0);
    nPeaks = 0;
    nPeaksArcFM = 0;
    nPeaksArcSep = 0;
    doPeakDetection = 0;
    peakMethod = "skipped";
    detectionType = "peaks";
    isRegionMethodUsed = 0;
    regionMinAreaUm = def_regionMinAreaUm;
    regionMaxAreaUm = def_regionMaxAreaUm;
    regionSmoothUm = def_regionSmoothUm;
    regionSepMode = "Intensity watershed (maxima)";
    regionPromPct = def_regionPromPct;
    regionDogUm = def_regionDogUm;
    peakPreviewZoomPct = 100;
    profileArcExport = newArray(0);
    arcValidExport = newArray(0);
    swArcExport = 0;
    arcProfMin = 0;
    arcProfMax = 0;
    arcProfValidCols = 0;
    peaksBright = true;
    tolPct = def_tolPct;
    minSepUm = def_minSepUm;
    markerUm = def_markerUm;
    doRefine = def_doRefine;
    requireWindow(analysisCorrWin);
    getDimensions(swPk, shPk, cPk, sPk, fPk);

    adjacentExportsDone = 0;
    masterAnalysisChIdx = 0;
    if (adjacentMode) {
        // Shared mask/geometry already built; run one tuning/export pass per analysis channel.
        runAdjacentDetectionPasses();
    } else if (batchMode) {
        if (def_usePeakDetection != 0) {
            doPeakDetection = 1;
            runAutoPeakDetection();
            requireWindow(maskTitle);
            beginHiddenCompute();
            finalizePeakDetectionResults();
            endHiddenCompute();
        } else {
            doPeakDetection = 0;
            peakMethod = "skipped";
            detectionType = "skipped";
            closeWindowsWithTitlePrefix(peakPreviewPrefix);
            closeWindowIfOpen(peakProfilePlotWin);
        }
    } else {
    peakDetectDefault = "No - skip";
    if (def_usePeakDetection != 0) peakDetectDefault = "Yes - tune and detect";
    Dialog.create("Patch detection");
    Dialog.addMessage("Core measurements are complete. Optionally detect patches / blobs along the neurite.");
    Dialog.addChoice("Perform patch detection?", newArray("Yes - tune and detect", "No - skip"), peakDetectDefault);
    Dialog.show();
    peakChoice = Dialog.getChoice();

    if (indexOf(peakChoice, "Yes") >= 0) {
        doPeakDetection = 1;
        prefSet("unified.doPeakDetection", 1);
        peakMethod = "Masked regions (threshold + particles)";
        prefSet("unified.peakMethod", peakMethod);
        requireWindow(maskTitle); getDimensions(swMpk, shMpk, cMpk, sMpk, fMpk);
        swArc = swPk;
        if (swMpk < swArc) swArc = swMpk;
        shArc = shPk;
        if (shMpk < shArc) shArc = shMpk;
        profileArc = newArray(swArc);
        yAtMax = newArray(swArc);
        arcValid = newArray(swArc);
        peakBufSize = swArc * shArc;
        if (peakBufSize < 1) peakBufSize = 1;
        xsWork = newArray(peakBufSize);
        ysWork = newArray(peakBufSize);
        valsWork = newArray(peakBufSize);
        peakRegAreaWork = newArray(peakBufSize);
        peakRegMeanWork = newArray(peakBufSize);
        peakRegBXWork = newArray(peakBufSize);
        peakRegBYWork = newArray(peakBufSize);
        peakRegWWork = newArray(peakBufSize);
        peakRegHWork = newArray(peakBufSize);
        peakRegFeretWork = newArray(peakBufSize);
        peakRegMinFeretWork = newArray(peakBufSize);
        peakRegCircWork = newArray(peakBufSize);
        peakRegSrcIdxWork = newArray(peakBufSize);
        peakRegDistMaskWork = newArray(peakBufSize);
        peakRegDistMaskCentWork = newArray(peakBufSize);
        peakRegDistMidWork = newArray(peakBufSize);
        peakCounts = newArray(3);
        peakRegionKeptRoiIdx = newArray(peakBufSize);
        regionDiagBXWork = newArray(peakBufSize);
        regionDiagBYWork = newArray(peakBufSize);
        regionDiagWWork = newArray(peakBufSize);
        regionDiagHWork = newArray(peakBufSize);
        setBatchMode(false);
        isRegions = 1;
        peakTuneTitle = "Patch detection tuning";
        donePeaks = false;
        peakPreviewReady = 0;
        while (!donePeaks) {
                Dialog.createNonBlocking(peakTuneTitle);
                    Dialog.addNumber("Tolerance (% of range):", tolPct, 2, 5, "Min patch brightness; lower = more patches");
                    Dialog.addNumber("Min inter-peak distance (um):", minSepUm, 2, 5, "Drop peaks closer than this");
                    Dialog.addNumber("Min region area (um^2):", regionMinAreaUm, 3, 5, "Min patch size; lower = smaller blobs");
                    Dialog.addNumber("Max region area (um^2, 0=none):", regionMaxAreaUm, 1, 5, "0 = no maximum patch size");
                    Dialog.addNumber("Pre-threshold Gaussian sigma (um, 0=off):", regionSmoothUm, 2, 5, "Blur analysis channel before detect; lower = sharper");
                    promLabel = "Maxima prominence (% of in-mask range):";
                    Dialog.addNumber(promLabel, regionPromPct, 1, 5, "Split merged patches; lower = more splits");
                    flattenHintUnits = "Regions method only";
                    if (regionDogUm <= 0) flattenHintUnits = "Off - unlikely to split patches";
                    else flattenHintUnits = "Subtract neurite glow (~3-4 um); 0 = off";
                    Dialog.addNumber("Baseline flatten diameter (um, 0=off):", regionDogUm, 2, 5, flattenHintUnits);
                if (peakPreviewReady == 1) {
                        if (nPeaks > 0) {
                            totalRegArea = sumRegionAreas(nPeaks);
                            meanRegArea = totalRegArea / nPeaks;
                            countMsg = "Blobs: " + d2s(nPeaks, 0) +
                                "  |  Mean area: " + d2s(meanRegArea, 2) + " um^2" +
                                "  |  Total area: " + d2s(totalRegArea, 2) + " um^2";
                        } else {
                            countMsg = "Blobs: 0";
                        }
                    Dialog.addMessage(countMsg);
                }
                if (peakPreviewReady == 1) peakActionDefault = "Accept patches and continue";
                else peakActionDefault = "Update preview";
                Dialog.addChoice("Action:", newArray("Update preview", "Accept patches and continue"), peakActionDefault);
                Dialog.show();
                tolPct = Dialog.getNumber();
                minSepUm = Dialog.getNumber();
                    regionMinAreaUm = Dialog.getNumber();
                    regionMaxAreaUm = Dialog.getNumber();
                    regionSmoothUm = Dialog.getNumber();
                    regionPromPct = Dialog.getNumber();
                    regionDogUm = Dialog.getNumber();
                actionPeaks = Dialog.getChoice();
                prefSet("unified.tolPct", tolPct);
                prefSet("unified.minSepUm", minSepUm);
                    prefSet("unified.regionMinAreaUm", regionMinAreaUm);
                    prefSet("unified.regionMaxAreaUm", regionMaxAreaUm);
                    prefSet("unified.regionSmoothUm", regionSmoothUm);
                    prefSet("unified.regionPromPct", regionPromPct);
                    prefSet("unified.regionDogUm", regionDogUm);
                doPreviewUpdate = 0;
                doAcceptPeaks = 0;
                if (indexOf(actionPeaks, "Update") >= 0) doPreviewUpdate = 1;
                if (indexOf(actionPeaks, "Accept") >= 0) {
                    doAcceptPeaks = 1;
                    if (peakPreviewReady == 0) {
                        showStatus("Running detection before accept...");
                        doPreviewUpdate = 1;
                    }
                }
                if (doPreviewUpdate) {
                    requireWindow(maskTitle);
                    closeWindowsWithTitlePrefix(peakPreviewPrefix);
                    closeWindowIfOpen(peakProfilePlotWin);
                    roiManager("reset");
                    beginHiddenCompute();
                    runPeakDetection(peakMethod, maskTitle, analysisCorrWin, swArc, shArc, tolPct, minSepUm, peaksBright, doRefine, umPerPx, regionMinAreaUm, regionMaxAreaUm, regionSmoothUm, profileArc, yAtMax, arcValid, xsWork, ysWork, valsWork, peakCounts);
                    endHiddenCompute();
                    nPeaksArcFM = peakCounts[0];
                    nPeaksArcSep = peakCounts[1];
                    nPeaks = peakCounts[2];
                    xs = newArray(nPeaks);
                    ys = newArray(nPeaks);
                    vals = newArray(nPeaks);
                    for (i = 0; i < nPeaks; i++) {
                        xs[i] = xsWork[i];
                        ys[i] = ysWork[i];
                        vals[i] = valsWork[i];
                    }
                        if (nPeaks > 0) {
                            buildRegionPreview(analysisCorrWin, nPeaks, xs, ys, peakRegBXWork, peakRegBYWork, peakRegWWork, peakRegHWork, peakPreviewPrefix);
                        } else if (regionDiagN > 0) {
                            buildRegionDiagnosticPreview(analysisCorrWin, regionDiagN, regionDiagBXWork, regionDiagBYWork, regionDiagWWork, regionDiagHWork, peakPreviewPrefix);
                        } else {
                            buildRegionPreview(analysisCorrWin, 0, xs, ys, peakRegBXWork, peakRegBYWork, peakRegWWork, peakRegHWork, peakPreviewPrefix);
                        }
                    zoomPeakPreviewForTuning(peakPreviewPrefix);
                    layoutPeakPreviews();
                    selectWindowByTitlePrefix(peakPreviewPrefix);
                    peakPreviewReady = 1;
                    if (nPeaks == 0 && regionDiagN > 0) {
                        showStatus("Diagnostic preview (orange): " + d2s(regionDiagN, 0) + " speckles below min area — raise Tolerance or lower Min region area.");
                    }
                }
                if (doAcceptPeaks) {
                    closeWindowsWithTitlePrefix(peakPreviewPrefix);
                    closeWindowIfOpen(peakProfilePlotWin);
                        isRegionMethodUsed = 1;
                        detectionType = "regions";
                        regAreaExport = newArray(nPeaks);
                        regMeanExport = newArray(nPeaks);
                        regBXExport = newArray(nPeaks);
                        regBYExport = newArray(nPeaks);
                        regWExport = newArray(nPeaks);
                        regHExport = newArray(nPeaks);
                        regFeretExport = newArray(nPeaks);
                        regMinFeretExport = newArray(nPeaks);
                        regCircExport = newArray(nPeaks);
                        regDistMaskExport = newArray(nPeaks);
                        regDistMaskCentExport = newArray(nPeaks);
                        regDistMidExport = newArray(nPeaks);
                        for (i = 0; i < nPeaks; i++) {
                            regAreaExport[i] = peakRegAreaWork[i];
                            regMeanExport[i] = peakRegMeanWork[i];
                            regBXExport[i] = peakRegBXWork[i];
                            regBYExport[i] = peakRegBYWork[i];
                            regWExport[i] = peakRegWWork[i];
                            regHExport[i] = peakRegHWork[i];
                            regFeretExport[i] = peakRegFeretWork[i];
                            regMinFeretExport[i] = peakRegMinFeretWork[i];
                            regCircExport[i] = peakRegCircWork[i];
                            regDistMaskExport[i] = peakRegDistMaskWork[i];
                            regDistMaskCentExport[i] = peakRegDistMaskCentWork[i];
                            regDistMidExport[i] = peakRegDistMidWork[i];
                        }
                    donePeaks = true;
                }
        }
        requireWindow(maskTitle);
        beginHiddenCompute();
        finalizePeakDetectionResults();
        endHiddenCompute();
        if (nPeaks == 0 && !batchMode) {
                showMessage("Patch detection", "No regions found at the current settings.\n" +
                    "Try a lower tolerance, smaller min area, or adjust the mask.\n" +
                    "After particles: " + d2s(nPeaksArcFM, 0) + " region(s).");
        }
    } else {
        closeWindowsWithTitlePrefix(peakPreviewPrefix);
        closeWindowIfOpen(peakProfilePlotWin);
        detectionType = "skipped";
        prefSet("unified.doPeakDetection", 0);
    }
    }

    // Save detection table (patches for masked regions; peaks for Find Maxima methods)
    // Adjacent mode already wrote per-channel exports inside runAdjacentDetectionPasses.
    if (doPeakDetection && !adjacentExportsDone) {
        if (isRegionMethodUsed) {
            savePatchesTable(nPeaks, regAreaExport, regMeanExport, regCircExport, xs, ys, regBXExport, regBYExport, regWExport, regHExport, regFeretExport, regMinFeretExport, swPk, roiLenUm, umPerPx, analysisOutDir + File.separator + base + analysisExportTag + "_patches.csv");
        } else {
            Table.create("Peaks");
            for (i = 0; i < nPeaks; i++) {
                // Arc length coordinate
                if (swPk > 1) xNormPk = xs[i] / (swPk - 1); else xNormPk = 0;
                xumPk = xNormPk * roiLenUm;
                Table.set("X px straight", i, xs[i]);
                Table.set("Y px straight", i, ys[i]);
                Table.set("X arc (um)", i, xumPk);
                Table.set("Peak value (AU)", i, vals[i]);
            }
            Table.save(analysisOutDir + File.separator + base + analysisExportTag + "_peaks.csv");
            if (isOpen("Peaks")) { selectWindow("Peaks"); run("Close"); }
            if (nPeaks > 1) {
                Table.create("PeakDistances");
                if (swPk > 1) px_to_um = roiLenUm / (swPk - 1); else px_to_um = 0;
                sortedX = sortValuesAscending(xs, nPeaks);
                for (i = 0; i < nPeaks - 1; i++) {
                    dpx = sortedX[i+1] - sortedX[i];
                    arcDist = dpx * px_to_um;
                    Table.set("Start X (um)", i, sortedX[i] * px_to_um);
                    Table.set("End X (um)", i, sortedX[i+1] * px_to_um);
                    Table.set("Distance (um)", i, arcDist);
                }
                Table.save(analysisOutDir + File.separator + base + analysisExportTag + "_peak_distances.csv");
                if (isOpen("PeakDistances")) { selectWindow("PeakDistances"); run("Close"); }
            }
        }
    }
    // Peak density per um^2: compute mask area
    totalMaskArea = 0;
    for (x = 0; x < widthsUm.length; x++) {
        // Each column corresponds to length along arc of umPerPx
        totalMaskArea += widthsUm[x] * umPerPx;
    }
    if (totalMaskArea > 0) peakDensity = nPeaks / totalMaskArea; else peakDensity = 0;
    // Particles adjacent to the neurite are conventionally reported per unit of neurite length,
    // since the mask area they are compared against does not contain them.
    if (roiLenUm > 0) adjLinearDensity = nPeaks / roiLenUm; else adjLinearDensity = 0;

    // ----- Save profiles -----
    if (includeCh2Profile) {
        minLen = ch2_inMask_raw.length;
        if (ch2_inMask_bgsub.length < minLen) minLen = ch2_inMask_bgsub.length;
        if (minLen < 0) minLen = 0;
        Table.create("Analysis Profile");
        for (i = 0; i < minLen; i++) {
            if (swC2c > 1) xnormPr = i / (swC2c - 1); else xnormPr = 0;
            xumPr   = xnormPr * roiLenUm;
            Table.set("X arc (um)", i, xumPr);
            Table.set("X norm 0-1", i, xnormPr);
            Table.set("Analysis Channel Intensity Raw (AU)", i, ch2_inMask_raw[i]);
            Table.set("Analysis Channel Intensity Bgsub (AU)", i, ch2_inMask_bgsub[i]);
            if (doPeakDetection && isArcPeakMethod(peakMethod)) {
                if (i < swArcExport) {
                    Table.set("Arc Valid", i, arcValidExport[i]);
                    if (arcValidExport[i]) {
                        Table.set("Arc Max Analysis Channel Bgsub (AU)", i, profileArcExport[i]);
                    } else {
                        Table.set("Arc Max Analysis Channel Bgsub (AU)", i, "");
                    }
                } else {
                    Table.set("Arc Valid", i, 0);
                    Table.set("Arc Max Analysis Channel Bgsub (AU)", i, "");
                }
                pkAt = 0;
                for (p = 0; p < nPeaks; p++) {
                    if (round(xs[p]) == i) pkAt = 1;
                }
                Table.set("Peak At Column", i, pkAt);
            }
        }
        Table.save(analysisOutDir + File.separator + base + analysisExportTag + "_profile.csv");
        if (isOpen("Analysis Profile")) { selectWindow("Analysis Profile"); run("Close"); }
    }

    if (includeRoiProfile) {
        minLenC1 = profileC1Raw.length;
        if (profileC1Bg.length < minLenC1) minLenC1 = profileC1Bg.length;
        Table.create("ROI Profile");
        for (i = 0; i < minLenC1; i++) {
            if (swC1 > 1) xnormC1 = i / (swC1 - 1); else xnormC1 = 0;
            xumC1   = xnormC1 * roiLenUm;
            Table.set("X px straight", i, i);
            Table.set("X arc (um)", i, xumC1);
            Table.set("X norm 0-1", i, xnormC1);
            Table.set("ROI Channel Intensity Raw (AU)", i, profileC1Raw[i]);
            Table.set("ROI Channel Intensity Bgsub (AU)", i, profileC1Bg[i]);
        }
        Table.save(roiOutDir + File.separator + base + roiExportTag + "_profile.csv");
        if (isOpen("ROI Profile")) { selectWindow("ROI Profile"); run("Close"); }
    }

    // ----- Summary tables -----
    if (includeRoiSummary) {
        Table.create("ROI Summary");
        Table.set("ROI Channel Length (um)", 0, roiLenUm);
        Table.set("ROI Channel Mean Intensity Raw (AU)", 0, meanC1Raw);
        Table.set("ROI Channel Mean Intensity Bgsub (AU)", 0, meanC1Bg);
        Table.set("ROI Channel Projection Integrated Intensity Raw (AU)", 0, projIntC1Raw);
        Table.set("ROI Channel Projection Integrated Intensity Bgsub (AU)", 0, projIntC1Bg);
        Table.set("ROI Channel Intensity CV", 0, roiCvGlobal);
        Table.set("Neurite Width CV", 0, maskWidthCv);
        if (includeMaskWidth) {
            Table.set("Neurite Width Mean (um)", 0, maskWidthMeanUm);
        }
        Table.save(roiOutDir + File.separator + base + roiExportTag + "_summary.csv");
        if (isOpen("ROI Summary")) { selectWindow("ROI Summary"); run("Close"); }
    }

    if (includeNP) {
        Table.create("Neuropil Summary");
        Table.set("Analysis Channel Mean Neuropil Raw (AU)", 0, ch2_mean_neuropil_raw);
        Table.set("Analysis Channel Mean Corrected Raw (AU)", 0, ch2_mean_corrected_raw);
        Table.set("Analysis Channel Mean Enrichment Raw", 0, ch2_mean_enrich_raw);
        Table.set("Analysis Channel Mean Neuropil Bgsub (AU)", 0, ch2_mean_neuropil_bgsub);
        Table.set("Analysis Channel Mean Corrected Bgsub (AU)", 0, ch2_mean_corrected_bgsub);
        Table.set("Analysis Channel Mean Enrichment Bgsub", 0, ch2_mean_enrich_bgsub);
        Table.save(analysisOutDir + File.separator + base + analysisExportTag + "_neuropil_summary.csv");
        if (isOpen("Neuropil Summary")) { selectWindow("Neuropil Summary"); run("Close"); }
    }

    if (!adjacentExportsDone && (includeCh2IntSummary || includeCvSummary || doPeakDetection)) {
        Table.create("Analysis Summary");
        if (includeCh2IntSummary) {
            Table.set("Analysis Channel Mean Intensity Raw (AU)", 0, ch2_mean_inMask_raw);
            Table.set("Analysis Channel Mean Intensity Bgsub (AU)", 0, ch2_mean_inMask_bgsub);
            Table.set("Analysis Channel Projection Integrated Intensity Raw (AU)", 0, projIntC2Raw);
            Table.set("Analysis Channel Projection Integrated Intensity Bgsub (AU)", 0, projIntC2Bg);
        }
        if (includeCvSummary) {
            Table.set("Analysis Channel Intensity CV", 0, cv_global);
            if (doDetrend) {
                Table.set("Analysis Channel Intensity CV Detrended", 0, cv_detrended);
            }
        }
        if (doPeakDetection) {
            if (isRegionMethodUsed) {
                Table.set("Analysis Patch Count", 0, nPeaks);
                Table.set("Analysis Patch Density per (um2)", 0, peakDensity);
                Table.set("Analysis Patch Mean Area (um2)", 0, meanOfArray(regAreaExport, nPeaks));
                Table.set("Analysis Patch Mean Intensity (AU)", 0, meanOfArray(regMeanExport, nPeaks));
                Table.set("Analysis Patch Maximum Diameter Mean (um)", 0, meanOfArray(regFeretExport, nPeaks));
                Table.set("Analysis Patch Minimum Diameter Mean (um)", 0, meanOfArray(regMinFeretExport, nPeaks));
                Table.set("Analysis Patch Circularity Mean", 0, meanOfArray(regCircExport, nPeaks));
                Table.set("Analysis Patch Equivalent Circle Diameter Mean (um)", 0, meanPatchEquivDiamUm(regAreaExport, nPeaks));
                if (adjacentMode) {
                    Table.set("Analysis Patch Linear Density (patches / um)", 0, adjLinearDensity);
                    Table.set("Analysis Patch Mean Distance From Patch Edge To Mask (um)", 0, meanOfArray(regDistMaskExport, nPeaks));
                    Table.set("Analysis Patch Mean Distance From Patch Centre To Midline (um)", 0, meanAbsOfArray(regDistMidExport, nPeaks));
                }
            } else {
                Table.set("Analysis Peak Count", 0, nPeaks);
                Table.set("Analysis Peak Density per (um2)", 0, peakDensity);
            }
        }
        Table.save(analysisOutDir + File.separator + base + analysisExportTag + "_summary.csv");
        if (isOpen("Analysis Summary")) { selectWindow("Analysis Summary"); run("Close"); }
    }

    if (doPeakDetection && isArcPeakMethod(peakMethod)) {
        Table.create("Arc Profile Summary");
        Table.set("Arc Profile Min Analysis Channel Bgsub (AU)", 0, arcProfMin);
        Table.set("Arc Profile Max Analysis Channel Bgsub (AU)", 0, arcProfMax);
        Table.set("Arc Profile Valid Columns", 0, arcProfValidCols);
        Table.save(analysisOutDir + File.separator + base + analysisExportTag + "_arc_profile_summary.csv");
        if (isOpen("Arc Profile Summary")) { selectWindow("Arc Profile Summary"); run("Close"); }
    }

    if (!adjacentExportsDone && appendMasterSummary && masterSummaryPath != "") {
        masterRow = buildMasterSummaryRow();
        appendMasterSummaryRow(masterSummaryPath, masterRow);
    }

    // Analysis parameters summary
    if (includeAnalysisParams) {
        Table.create("Analysis Parameters");
        Table.set("ROI Channel Length (um)", 0, roiLenUm);
        Table.set("ROI Channel Index", 0, roiIdx);
        if (adjacentMode) {
            Table.set("Analysis Ch1 Channel Index", 0, ch2Idx);
        } else {
            Table.set("Analysis Channel Index", 0, ch2Idx);
        }
        Table.set("Mask Channel Index", 0, maskIdx);
        Table.set("FOV Width (um)", 0, fovUm);
        Table.set("Straighten Width (um)", 0, straightenUm);
        Table.set("Rolling Ball Radius (um)", 0, rbUm);
        Table.set("Mask Threshold Fraction", 0, maskThrFrac);
        Table.set("Mask Dilation (um)", 0, maskDilUm);
        Table.set("Tissue Type", 0, tissueType);
        if (adjacentMode) {
            Table.set("Adjacent Particle Mode", 0, adjacentMode);
            if (adjDualMode) Table.set("Analysis Ch2 Channel Index", 0, adjSecondAnalysisIdx);
        }
        if (is3DProject) {
            Table.set("Z Projection First Slice", 0, projZFirst);
            Table.set("Z Projection Last Slice", 0, projZLast);
            Table.set("Z Projection Method", 0, zProjMethod);
        }
        if (is3DPerPlane) {
            Table.set("3D Segment Count", 0, pathSegCount);
            Table.set("Compensate Z Segment Intensity", 0, compensateZIntensity);
        }
        Table.set("Include Analysis Parameters Export", 0, includeAnalysisParams);
        Table.set("Include Analysis Channel Intensity Summary", 0, includeCh2IntSummary);
        Table.set("Include Analysis Channel Profile", 0, includeCh2Profile);
        Table.set("Include Analysis Channel Intensity CV", 0, includeCvSummary);
        if (includeCv) {
            Table.set("Detrended CV", 0, doDetrend);
            if (doDetrend) Table.set("Detrend Window (um)", 0, detrendUm);
        }
        Table.set("Include Neurite Width", 0, includeMaskWidth);
        Table.set("Include ROI Channel Summary", 0, includeRoiSummary);
        Table.set("Include ROI Channel Profile", 0, includeRoiProfile);
        Table.set("Include ROI Channel Intensity", 0, includeCh1);
        Table.set("Include Neuropil", 0, includeNP);
        if (doPeakDetection) {
            Table.set("Peak Detection Performed", 0, doPeakDetection);
            if (adjacentMode) {
                Table.set("Analysis Ch1 Tolerance Percent", 0, adjParamTol1);
                Table.set("Analysis Ch1 Min Peak Separation (um)", 0, adjParamMinSep1);
                Table.set("Analysis Ch1 Region Min Area (um2)", 0, adjParamMinArea1);
                Table.set("Analysis Ch1 Region Max Area (um2)", 0, adjParamMaxArea1);
                Table.set("Analysis Ch1 Region Gaussian Sigma (um)", 0, adjParamSmooth1);
                Table.set("Analysis Ch1 Region Maxima Prominence Percent", 0, adjParamProm1);
                Table.set("Analysis Ch1 Region Baseline Flatten (um)", 0, adjParamDog1);
                Table.set("Adjacent Max Distance From Midline (um)", 0, adjMaxDistMidUm);
                if (adjDualMode) {
                    Table.set("Analysis Ch2 Tolerance Percent", 0, adjParamTol2);
                    Table.set("Analysis Ch2 Min Peak Separation (um)", 0, adjParamMinSep2);
                    Table.set("Analysis Ch2 Region Min Area (um2)", 0, adjParamMinArea2);
                    Table.set("Analysis Ch2 Region Max Area (um2)", 0, adjParamMaxArea2);
                    Table.set("Analysis Ch2 Region Gaussian Sigma (um)", 0, adjParamSmooth2);
                    Table.set("Analysis Ch2 Region Maxima Prominence Percent", 0, adjParamProm2);
                    Table.set("Analysis Ch2 Region Baseline Flatten (um)", 0, adjParamDog2);
                }
            } else {
                Table.set("Tolerance Percent", 0, tolPct);
                Table.set("Min Peak Separation (um)", 0, minSepUm);
                if (isArcPeakMethod(peakMethod)) {
                    Table.set("Peak Marker Diameter (um)", 0, markerUm);
                    Table.set("Refine Peaks", 0, doRefine);
                }
                if (isRegionMethodUsed) {
                    Table.set("Region Min Area (um2)", 0, regionMinAreaUm);
                    Table.set("Region Max Area (um2)", 0, regionMaxAreaUm);
                    Table.set("Region Gaussian Sigma (um)", 0, regionSmoothUm);
                    Table.set("Region Maxima Prominence Percent", 0, regionPromPct);
                    Table.set("Region Baseline Flatten (um)", 0, regionDogUm);
                }
            }
        }
        Table.save(paramsOutDir + File.separator + base + "_analysis_parameters.csv");
        if (isOpen("Analysis Parameters")) { selectWindow("Analysis Parameters"); run("Close"); }
    }

    // Plain straightened PNGs saved before peak tuning (saveStraightenedPlainPngs).

    // ----- Overlay peaks drawing -----
    if (doPeakDetection && !adjacentExportsDone) {
        beginHiddenCompute();
        activateImage(analysisCorrWin);
        duplicateCurrentAs(analysisOverlayWin);
        run("Remove Overlay"); setOption("Display Label", false);
        addPeakOverlaysToImage(isRegionMethodUsed, nPeaks, xs, ys, markerUm, umPerPx);
        run("Flatten");
        overlayPeaksPath = analysisOutDir + File.separator + base + analysisExportTag + "_overlay_patches.png";
        saveAs("PNG", overlayPeaksPath);
        closeWindowAfterSave(analysisOverlayWin, overlayPeaksPath);
        roiManager("reset");
        endHiddenCompute();
    }

    // ----- Cleanup windows -----
    // Close channel windows and intermediates
    titlesToClose = newArray(roiWin, analysisWin, analysis2Win, maskChannelWin, roiStraightWin, analysisStraightWin, analysis2StraightWin, roiBgsubPreviewWin, analysisBgsubPreviewWin, analysis2BgsubPreviewWin, "Mask_straight", roiCorrWin, analysisCorrWin, analysis2CorrWin, analysisNormWin, "ArcProfile_FM", maskPreviewWin, analysisMaskPreviewWin, analysis2MaskPreviewWin, peakProfilePlotWin, "Regions_proc", "Regions_binary", "Regions_binary_pre", "Regions_intensity", "Regions_measure", "Regions_src", "Regions_small", "Regions_base", "Mask_work", "Mask_seed", "Mask_final", "Mask_core", roiPlainWin, analysisPlainWin, analysis2PlainWin, analysisOverlayWin, analysis2OverlayWin, adjMergeRoiWin, adjMergedOverlayWin, maskTitle);
    for (i = 0; i < titlesToClose.length; i++) {
        t = titlesToClose[i];
        if (t != "" && isOpen(t)) { selectWindow(t); run("Close"); }
    }
    closeWindowsWithTitlePrefix(peakPreviewPrefix);
    if (adjDualMode && analysis2Win != "") closeWindowsWithTitlePrefix(analysis2Win + "_peak_preview");
    closeWindowsWithTitlePrefix("Result of ");
    closeWindowsWithTitlePrefix("Plot of ");
    closeWindowAfterSave("", roiOutDir + File.separator + base + roiExportTag + "_straight.png");
    closeWindowAfterSave("", analysisOutDir + File.separator + base + analysisExportTag + "_straight.png");
    if (adjDualMode && analysis2ExportTag != "") closeWindowAfterSave("", analysisOutDir + File.separator + base + analysis2ExportTag + "_straight.png");
    if (doPeakDetection) closeWindowAfterSave("", overlayPeaksPath);
    if (adjDualMode) {
        closeWindowAfterSave("", analysisOutDir + File.separator + base + "_overlay_patches_merged.png");
        closeWindowsWithTitlePrefix(adjMergedOverlayWin);
    }
    // Close any residual 'C1-' or 'C2-' split channel windows
    titlesAll = getList("image.titles");
    if (titlesAll.length > 0) {
        for (i = 0; i < titlesAll.length; i++) {
            t = titlesAll[i];
            if (startsWith(t, "C1-") || startsWith(t, "C2-") || indexOf(t, "Channel 1") != -1 || indexOf(t, "Channel 2") != -1) {
                if (isOpen(t)) { selectWindow(t); run("Close"); }
            }
        }
    }
    // Close any results windows
    if (isOpen("Results")) { selectWindow("Results"); run("Close"); }
    if (isOpen("ROI Manager")) { selectWindow("ROI Manager"); run("Close"); }
    batchProcessedCount++;
    if (!batchMode) {
        showMessage("Done", "Analysis complete. Results saved to:\n" + outDir);
        batchMore = false;
    } else {
        def_lastInputDir = File.getDirectory(inPath);
        if (File.exists(inPath) && !File.isDirectory(inPath)) def_batchInputSeed = inPath; else def_batchInputSeed = def_lastInputDir;
        Dialog.create("Neurite FluoroQuant - Fast Compute - Batch (next image)");
        Dialog.addMessage("Select the next image, or finish the batch.");
        Dialog.addFile("Open next image:", def_batchInputSeed);
        Dialog.addCheckbox("Finish batch", false);
        Dialog.show();
        inPath = Dialog.getString();
        finishBatch = Dialog.getCheckbox();
        if (finishBatch || inPath == "") batchMore = false;
    }
    }
    if (batchMode && batchProcessedCount > 0) {
        showMessage("Batch complete", "Processed " + d2s(batchProcessedCount, 0) + " image(s).\nResults in: " + outParent);
    }
}