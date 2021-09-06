if (app.documents.length > 0) {
    main();
} else {
    Window.alert("No active document found.");
}

function getNewFile(destFolder, filename, extension) {
    var newName = filename + "." + extension;
    var newFile = new File(destFolder + '/' + newName);
    return newFile;
}

function getIllustratorOptions() {
    var saveOptions = new IllustratorSaveOptions();
    saveOptions.compressed = true;
    saveOptions.embedICCProfile = true;
    saveOptions.embedLinkedFiles = true;
    saveOptions.pdfCompatible = true;
    saveOptions.saveMultipleArtboards = true;

    return saveOptions;
}

function getPDFOptions() {
    var saveOptions = new PDFSaveOptions();
    saveOptions.optimization = true;
    saveOptions.preserveEditability = false;
    saveOptions.embedICCProfile = true;
    saveOptions.viewAfterSaving = false;

    return saveOptions;
}

function getPNGOptions() {
    var saveOptions = new ExportOptionsPNG24();
    saveOptions.antiAliasingMethod = AntiAliasingMethod.TYPEOPTIMIZED;
    saveOptions.antiAliasing = true;
    saveOptions.transparency = true;
    saveOptions.artBoardClipping = false;
    saveOptions.matte = true;

    // resolution is not supported.
    // saveOptions.resolution = 300;
    saveOptions.verticalScale = 416.7;
    saveOptions.horizontalScale = 416.7;

    return saveOptions;
}

function getImageCaptureOptions() {
    var saveOptions = new ImageCaptureOptions();
    saveOptions.antiAliasing = true;
    saveOptions.antiAliasingMethod = AntiAliasingMethod.TYPEOPTIMIZED;
    saveOptions.artBoardClipping = false;
    saveOptions.matte = true;
    saveOptions.resolution = 300;
    saveOptions.transparency = true;

    return saveOptions;
}

function openFile(file) {
    app.open(file);
    return app.activeDocument;
}
function main() {
    var document = app.activeDocument;
    var currentFile = document.fullName;
    var filename = document.name.split('.')[0];
    var destFolder = app.activeDocument.fullName.parent;
    var tmpFolder = new Folder("/tmp")

    var embedIllustratorFile = getNewFile(tmpFolder, filename, "ai");
    var aiSaveOptions = getIllustratorOptions();

    document.saveAs(embedIllustratorFile, aiSaveOptions);
    document.close(SaveOptions.DONOTSAVECHANGES);

    document = openFile(embedIllustratorFile);
    var pdfFile = getNewFile(destFolder, filename, "pdf");
    var pdfSaveOptions = getPDFOptions();
    document.saveAs(pdfFile, pdfSaveOptions);
    document.close(SaveOptions.DONOTSAVECHANGES);

    // document = openFile(embedIllustratorFile);
    // var pngFile = getNewFile(destFolder, filename, "png");
    // // var pngSaveOptions = getPNGOptions();
    // // document.exportFile(pngFile, ExportType.PNG24, pngSaveOptions);

    // var idx = document.artboards.getActiveArtboardIndex();
    // var activeArtboard = document.artboards[idx];

    // var imageCaptureOptions = getImageCaptureOptions();
    // document.imageCapture(pngFile, activeArtboard.artboardRect, imageCaptureOptions)
    // document.close(SaveOptions.DONOTSAVECHANGES);

    Window.alert("Exported to " + destFolder.fsName);

    app.open(currentFile);
}

