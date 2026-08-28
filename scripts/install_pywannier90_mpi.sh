#!/usr/bin/env bash
# Install MPI-enabled Wannier90-DLWF and its MPI-aware PySCF adapter.
#
# The pyWannier90 adapter remains BSD-3-Clause licensed and retains the
# original attribution in pyWannier90/LICENSE and THIRD_PARTY_NOTICES.md.

set -euo pipefail

usage() {
    cat <<'EOF'
Usage: scripts/install_pywannier90_mpi.sh [options]

Build and install the MPI-enabled Wannier90 v4 executable and shared C
library, install the mpi4py-aware pywannier90_v4 adapter, generate an
activation file, and run the single-rank integration tests.

Options:
  --prefix PATH       Installation prefix
                      [default: $HOME/.local/wannier90-dlwf-4.0-mpi]
  --build-dir PATH    CMake build directory
                      [default: REPOSITORY/build/build-pywannier90-mpi]
  --python COMMAND    Python containing PySCF 2.9+ and MPI-compatible mpi4py
                      [default: python]
  --jobs N            Parallel compilation jobs [default: detected]
  --skip-python       Do not install the Python adapter
  --skip-tests        Do not run the Python/native integration tests
  -h, --help          Show this help

MPI and compiler selection can be supplied through the usual environment or
CMake variables. For example:
  CC=gcc FC=mpifort scripts/install_pywannier90_mpi.sh

Load the same compiler and MPI module stack before building and running.

The resulting library is named libwannier90_mpi, not libwannier90. mpi4py must
be built against the same MPI implementation used to build Wannier90.
EOF
}

die() {
    echo "error: $*" >&2
    exit 1
}

require_value() {
    [[ $# -ge 2 ]] || die "$1 requires a value"
    [[ -n "${2:-}" ]] || die "$1 requires a value"
}

script_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)"
repository_root="$(cd -- "$script_dir/.." && pwd -P)"

install_prefix="${HOME}/.local/wannier90-dlwf-4.0-mpi"
build_dir="$repository_root/build/build-pywannier90-mpi"
python_command="python"
install_python=true
run_tests=true

if command -v nproc >/dev/null 2>&1; then
    jobs="$(nproc)"
elif command -v sysctl >/dev/null 2>&1; then
    jobs="$(sysctl -n hw.ncpu 2>/dev/null || echo 4)"
else
    jobs=4
fi

while [[ $# -gt 0 ]]; do
    case "$1" in
        --prefix)
            require_value "$@"
            install_prefix="$2"
            shift 2
            ;;
        --build-dir)
            require_value "$@"
            build_dir="$2"
            shift 2
            ;;
        --python)
            require_value "$@"
            python_command="$2"
            shift 2
            ;;
        --jobs)
            require_value "$@"
            jobs="$2"
            [[ "$jobs" =~ ^[1-9][0-9]*$ ]] || die "--jobs must be a positive integer"
            shift 2
            ;;
        --skip-python)
            install_python=false
            shift
            ;;
        --skip-tests)
            run_tests=false
            shift
            ;;
        -h|--help)
            usage
            exit 0
            ;;
        *)
            die "unknown option: $1 (use --help)"
            ;;
    esac
done

command -v cmake >/dev/null 2>&1 || die "cmake was not found in PATH"
command -v mpiexec >/dev/null 2>&1 || \
    die "mpiexec was not found; load an MPI implementation before installing"

if $install_python || $run_tests; then
    command -v "$python_command" >/dev/null 2>&1 || \
        die "Python interpreter was not found: $python_command"
    if [[ ! -f "$repository_root/pyWannier90/pyproject.toml" ]]; then
        if git -C "$repository_root" rev-parse --is-inside-work-tree >/dev/null 2>&1; then
            echo "Initializing the pyWannier90 submodule..."
            git -C "$repository_root" submodule update --init --recursive -- pyWannier90
        fi
    fi
    [[ -f "$repository_root/pyWannier90/pyproject.toml" ]] || \
        die "pyWannier90 is unavailable; run 'git submodule update --init --recursive'"
fi

mkdir -p -- "$build_dir" "$install_prefix"
build_dir="$(cd -- "$build_dir" && pwd -P)"
install_prefix="$(cd -- "$install_prefix" && pwd -P)"

echo "Configuring MPI-enabled Wannier90:"
echo "  source:  $repository_root"
echo "  build:   $build_dir"
echo "  prefix:  $install_prefix"
echo "  Python:  $python_command"
echo "  jobs:    $jobs"

cmake \
    -S "$repository_root" \
    -B "$build_dir" \
    -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_INSTALL_PREFIX="$install_prefix" \
    -DWANNIER90_WITH_C=ON \
    -DWANNIER90_SHARED_LIBS=ON \
    -DWANNIER90_INSTALL=ON \
    -DWANNIER90_TEST=OFF \
    -DWANNIER90_MPI=ON

cmake --build "$build_dir" --parallel "$jobs"
cmake --install "$build_dir"

installed_library=""
for candidate in \
    "$install_prefix/lib/libwannier90_mpi.so" \
    "$install_prefix/lib64/libwannier90_mpi.so" \
    "$install_prefix/lib/libwannier90_mpi.dylib" \
    "$install_prefix/lib64/libwannier90_mpi.dylib"; do
    if [[ -e "$candidate" ]]; then
        installed_library="$candidate"
        break
    fi
done

if [[ -z "$installed_library" ]]; then
    installed_library="$(find -L "$install_prefix" -type f \
        \( -name 'libwannier90_mpi.so' -o -name 'libwannier90_mpi.dylib' \) \
        -print -quit)"
fi
[[ -n "$installed_library" ]] || \
    die "the installed libwannier90_mpi shared library could not be found"

installed_executable="$install_prefix/bin/wannier90.x"
[[ -x "$installed_executable" ]] || \
    die "the installed MPI wannier90.x executable could not be found"

library_dir="$(dirname -- "$installed_library")"
activation_file="$install_prefix/activate-wannier90-mpi.sh"

cat >"$activation_file" <<EOF
#!/usr/bin/env bash
# Generated by scripts/install_pywannier90_mpi.sh
# Load the same MPI/compiler module stack used for the build before sourcing.
export WANNIER90_MPI_ROOT="$install_prefix"
export WANNIER90_MPI_LIB="$installed_library"
export PYWANNIER90_LIB="$installed_library"
export PATH="$install_prefix/bin\${PATH:+:\$PATH}"
export LD_LIBRARY_PATH="$library_dir\${LD_LIBRARY_PATH:+:\$LD_LIBRARY_PATH}"
export DYLD_LIBRARY_PATH="$library_dir\${DYLD_LIBRARY_PATH:+:\$DYLD_LIBRARY_PATH}"
export PKG_CONFIG_PATH="$library_dir/pkgconfig\${PKG_CONFIG_PATH:+:\$PKG_CONFIG_PATH}"
export CMAKE_PREFIX_PATH="$install_prefix\${CMAKE_PREFIX_PATH:+:\$CMAKE_PREFIX_PATH}"
EOF
chmod 755 "$activation_file"

if $install_python || $run_tests; then
    "$python_command" - <<'PY'
import re

import numpy
import pyscf
try:
    import mpi4py
    from mpi4py import MPI
except ImportError as exc:
    raise SystemExit(
        "mpi4py is required and must use the same MPI implementation as Wannier90"
    ) from exc

numbers = tuple(int(item) for item in re.findall(r"\d+", pyscf.__version__)[:2])
if numbers < (2, 9):
    raise SystemExit(f"PySCF 2.9+ is required; found {pyscf.__version__}")
print(f"Using NumPy {numpy.__version__} and PySCF {pyscf.__version__}")
print(f"Using mpi4py {mpi4py.__version__}: {MPI.Get_library_version().splitlines()[0]}")
PY
fi

if $install_python; then
    "$python_command" -m pip install \
        --no-deps \
        --editable "$repository_root/pyWannier90"
fi

if $install_python || $run_tests; then
    env \
        PYTHONPATH="$repository_root/pyWannier90/src${PYTHONPATH:+:$PYTHONPATH}" \
        WANNIER90_MPI_LIB="$installed_library" \
        PYWANNIER90_LIB="$installed_library" \
        "$python_command" - <<'PY'
from pywannier90_v4 import Wannier90V4Library

with Wannier90V4Library() as library:
    library.require_site_symmetry_support()
    if not library.is_mpi:
        raise SystemExit(f"Expected an MPI Wannier90 library, got {library.path}")
    print(
        f"Loaded MPI Wannier90 v4 C library on rank "
        f"{library.mpi_rank}/{library.mpi_size}: {library.path}"
    )
PY
fi

if $run_tests; then
    echo "Running single-rank pyWannier90 MPI-library tests..."
    env \
        PYTHONPATH="$repository_root/pyWannier90/src${PYTHONPATH:+:$PYTHONPATH}" \
        WANNIER90_MPI_LIB="$installed_library" \
        PYWANNIER90_LIB="$installed_library" \
        PYWANNIER90_TEST_LIB="$installed_library" \
        "$python_command" -m unittest discover \
            -s "$repository_root/pyWannier90/tests" -v
fi

echo
echo "MPI-enabled Wannier90 installation complete."
echo "Wannier90 executable: $installed_executable"
echo "Wannier90 MPI library: $installed_library"
echo "Activate in the current shell with:"
echo "  source \"$activation_file\""
echo "Run a standalone calculation with, for example:"
echo "  mpiexec -n 4 \"$installed_executable\" seedname"
if $install_python || $run_tests; then
    echo "Run the bundled MPI PySCF/DLWF example with:"
    echo "  mpiexec -n 4 \"$python_command\" \"$repository_root/pyWannier90/examples/PySCF/h2_v4.py\""
fi
echo "Keep the MPI/compiler module stack used for the build loaded at runtime."
