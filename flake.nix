{
  description = "Nix dev shell for coord.serial";

  inputs.nixpkgs.url = "github:NixOS/nixpkgs/nixos-unstable";

  outputs = { self, nixpkgs }:
    let
      systems = [ "aarch64-darwin" "x86_64-darwin" "x86_64-linux" "aarch64-linux" ];
      forAllSystems = f:
        nixpkgs.lib.genAttrs systems (system:
          f (import nixpkgs { inherit system; }));
    in {
      devShells = forAllSystems (pkgs: {
        default = pkgs.mkShell {
          packages = with pkgs; [
            rWrapper
            rPackages.ggplot2
            rPackages.rlang
            rPackages.testthat
            rPackages.roxygen2
            rPackages.devtools
          ];

          shellHook = ''
            export R_LIBS_USER="''${XDG_DATA_HOME:-$HOME/.local/share}/coord_serial-nix/R-library"
            mkdir -p "$R_LIBS_USER"

            echo "Nix shell ready for coord.serial"
            echo "R_LIBS_USER=$R_LIBS_USER"
            echo "Run: ./nix/test-package.sh"
          '';
        };
      });
    };
}
