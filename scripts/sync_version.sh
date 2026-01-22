#!/bin/bash
# sync_version.sh - Synchronize version across project files
#
# This script updates the version number in:
# - VERSION file
# - Dockerfile (LABEL version and comment)
# - README.md (version references)
#
# Usage: ./scripts/sync_version.sh <version>
# Example: ./scripts/sync_version.sh 1.2.2

set -euo pipefail

VERSION="${1:-}"
if [ -z "$VERSION" ]; then
    echo "Error: Version argument required" >&2
    echo "Usage: $0 <version>" >&2
    exit 1
fi

# Remove 'v' prefix if present (semantic-release provides version without 'v')
VERSION_CLEAN="${VERSION#v}"

# Get script directory
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"

cd "$PROJECT_ROOT"

echo "Syncing version to ${VERSION_CLEAN}..."

# 1. Update VERSION file
echo "${VERSION_CLEAN}" > VERSION
echo "✓ Updated VERSION file"

# 2. Update Dockerfile
if [ -f "Dockerfile" ]; then
    # Update version comment
    sed -i "s/^# Version: .*/# Version: ${VERSION_CLEAN}/" Dockerfile
    
    # Update LABEL version
    sed -i "s/LABEL version=\".*\"/LABEL version=\"${VERSION_CLEAN}\"/" Dockerfile
    
    echo "✓ Updated Dockerfile"
fi

# 3. Update README.md
if [ -f "README.md" ]; then
    # Update version in **Version**: line
    sed -i "s/\*\*Version\*\*: [0-9.]*/**Version**: ${VERSION_CLEAN}/" README.md
    
    echo "✓ Updated README.md"
fi

echo "Version synchronization complete: ${VERSION_CLEAN}"
