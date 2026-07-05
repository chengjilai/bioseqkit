{
  description = "bioseqkit: a lightweight, dependency-free biological sequence processing toolkit";

  inputs = {
    nixpkgs.url = "git+ssh://git@github.com/NixOS/nixpkgs.git?ref=nixos-unstable&shallow=1";
    flake-utils.url = "git+ssh://git@github.com/numtide/flake-utils.git?shallow=1";
  };

  outputs =
    { self, nixpkgs, flake-utils }:
    flake-utils.lib.eachDefaultSystem (
      system:
      let
        pkgs = import nixpkgs { inherit system; };
        python = pkgs.python312;

        # The Python package itself, built from this repository via hatchling.
        # The core library has no runtime dependencies, so the closure is tiny.
        bioseqkit = python.pkgs.buildPythonPackage {
          pname = "bioseqkit";
          version = "0.2.0";
          pyproject = true;
          src = ./.;

          build-system = [ python.pkgs.hatchling ];

          nativeCheckInputs = [ python.pkgs.pytestCheckHook ];
          pythonImportsCheck = [ "bioseqkit" ];

          meta = with pkgs.lib; {
            description = "Lightweight, dependency-free biological sequence toolkit";
            homepage = "https://pypi.org/project/bioseqkit/";
            license = licenses.mit;
            mainProgram = "bioseqkit";
          };
        };

        # Reproducible OCI/Docker image built from the same derivation.
        #   nix build .#docker && docker load < result
        dockerImage = pkgs.dockerTools.buildLayeredImage {
          name = "bioseqkit";
          tag = "0.2.0";
          contents = [ bioseqkit ];
          config = {
            Entrypoint = [ "${bioseqkit}/bin/bioseqkit" ];
            Cmd = [ "--help" ];
          };
        };

        # Reproducible Apptainer/Singularity image (SIF) from the same derivation.
        #   nix build .#apptainer  (produces a .img/.sif)
        apptainerImage = pkgs.singularity-tools.buildImage {
          name = "bioseqkit-0.2.0";
          contents = [ bioseqkit ];
          runScript = "exec ${bioseqkit}/bin/bioseqkit \"$@\"";
          diskSize = 1024;
          memSize = 1024;
        };
      in
      {
        packages = {
          default = bioseqkit;
          bioseqkit = bioseqkit;
          docker = dockerImage;
          apptainer = apptainerImage;
        };

        # `nix run .` -> runs the bioseqkit CLI
        apps.default = flake-utils.lib.mkApp {
          drv = bioseqkit;
          name = "bioseqkit";
        };

        # `nix develop` -> fully reproducible development environment
        devShells.default = pkgs.mkShell {
          packages = [
            (python.withPackages (
              ps: with ps; [
                pytest
                pytest-benchmark
                ruff
                matplotlib
                seaborn
                jupyter
                sphinx
                myst-parser
              ]
            ))
            pkgs.uv
            pkgs.typst
          ];
        };

        checks.default = bioseqkit;
      }
    );
}
