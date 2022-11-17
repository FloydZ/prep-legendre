with import <nixpkgs> {};
let
  my-python = pkgs.python3;
  python-with-my-packages = my-python.withPackages (p: with p; [
  scipy
  ]);
in
{ pkgs ? import <nixpkgs> {} }:

stdenv.mkDerivation {
  name = "prep-legendre";
  src = ./.;

  buildInputs = [ python-with-my-packages git gnumake cmake clang ];

   shellHook = ''
	  # python hook
      PYTHONPATH=${python-with-my-packages}/${python-with-my-packages.sitePackages}
   ''; 
}
