{pkgs}: {
  deps = [
    pkgs.tk
    pkgs.tcl
    pkgs.qhull
    pkgs.gtk3
    pkgs.gobject-introspection
    pkgs.ghostscript
    pkgs.freetype
    pkgs.ffmpeg-full
    pkgs.cairo
    pkgs.xsimd
    pkgs.pkg-config
    pkgs.libxcrypt
  ];
  env = {
    PYTHON_LD_LIBRARY_PATH = pkgs.lib.makeLibraryPath [
      pkgs.tk
      pkgs.tcl
      pkgs.qhull
      pkgs.gtk3
      pkgs.gobject-introspection
      pkgs.ghostscript
      pkgs.freetype
      pkgs.cairo
      pkgs.xsimd
      pkgs.libxcrypt
    ];
  };
}
