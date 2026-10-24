#!/usr/bin/env bash
# SPDX-License-Identifier: Apache-2.0
#
# Stage premake for macOS builds. Packman publishes no macOS premake, so this puts premake's
# own 5.0.0-beta2 release (the version the Linux and Windows builds pin) into the Packman
# cache as the package deps/host-deps.packman.xml names. Its macOS binary is x86_64 and runs
# under Rosetta on Apple Silicon; it only writes the makefiles.
#
# Run once before ./repo.sh build -p macos-universal.
set -euo pipefail

VERSION="5.0.0-beta2"
URL="https://github.com/premake/premake-core/releases/download/v${VERSION}/premake-${VERSION}-macosx.tar.gz"
SHA256="620778e24847d4f8e2380cd98922977afec77cc8e805a77edcae1c05e0f6d44a"

if [[ "$(uname -s)" != "Darwin" ]]; then
    echo "setup-macos-premake.sh is for macOS only." >&2
    exit 1
fi

CACHE="${PM_PACKAGES_ROOT:-${HOME}/Library/Application Support/packman-cache}"
PACKAGE="${CACHE}/chk/premake/${VERSION}-macos-universal"

if [[ -x "${PACKAGE}/bin/premake5" ]]; then
    echo "premake ${VERSION} is already staged at ${PACKAGE}"
    exit 0
fi

WORK="$(mktemp -d)"
trap 'rm -rf "${WORK}"' EXIT

curl -sSLf -o "${WORK}/premake.tar.gz" "${URL}"
echo "${SHA256}  ${WORK}/premake.tar.gz" | shasum -a 256 -c - > /dev/null
tar -xzf "${WORK}/premake.tar.gz" -C "${WORK}" premake5

mkdir -p "${PACKAGE}/bin"
install -m 0755 "${WORK}/premake5" "${PACKAGE}/bin/premake5"
# Packman counts a package as installed when these are present.
touch "${PACKAGE}/.packman.lock" "${PACKAGE}/.packman.sha1"

echo "Staged premake ${VERSION} at ${PACKAGE}"
