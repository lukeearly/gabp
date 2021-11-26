with import <nixpkgs> {};

stdenv.mkDerivation {
  name = "gabp";
  version = "0.0.1";

  src = ../.;

  nativeBuildInputs = [ clang cmake ];

  cmakeFlags = [
    "-DBUILD_TESTING=1"
    #"-DCMAKE_BUILD_TYPE=Debug"
  ];

  buildPhase = "make -j $NIX_BUILD_CORES";

  installPhase = ''
    mkdir -p $out/bin
    #cp $TMP/gabp/build/test/gabp-tests $out/bin/
    cp -r $TMP $out
  '';
}
