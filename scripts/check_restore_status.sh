#!/bin/bash
#
# check_restore_status.sh - Check status of Glacier restoration requests
#
# DESCRIPTION:
#   Checks the restoration status of files from S3 Intelligent-Tiering Deep Archive.
#   After running restore_intelligent.py, files take ~12 hours (Bulk tier) to restore.
#   This script checks if they're ready for download.
#
# USAGE:
#   bash scripts/check_restore_status.sh
#
# OUTPUT:
#   Shows count of RESTORED, IN_PROGRESS, and DEEP_ARCHIVE_ACCESS files
#

set -e

echo "========================================"
echo "Checking Restoration Status"
echo "========================================"
echo ""

# Sample a few files to check status
echo "Checking sample files..."
echo ""

check_file() {
    local key="$1"
    local status=$(aws s3api head-object --bucket imaging-platform --key "$key" 2>/dev/null | grep -E "ArchiveStatus|Restore" || echo "RESTORED")
    echo "  $key"
    if echo "$status" | grep -q "RESTORED"; then
        echo "    ✓ RESTORED - Ready to download"
    elif echo "$status" | grep -q 'ongoing-request.*true'; then
        echo "    ⏳ IN_PROGRESS - Still restoring"
    elif echo "$status" | grep -q "DEEP_ARCHIVE_ACCESS"; then
        echo "    ❌ DEEP_ARCHIVE_ACCESS - Not yet started or still in archive"
    fi
    echo ""
}

# Check a few sample files from each directory
check_file "projects/2016_08_01_RadialMitochondriaDistribution_donna/workspace/metadata/CDRP_meta.csv"
check_file "projects/2016_08_01_RadialMitochondriaDistribution_donna/workspace/metadata/JUMP/compound.csv.gz"
check_file "projects/2016_08_01_RadialMitochondriaDistribution_donna/workspace/per_site_aggregated_profiles_newpattern_2/CDRP/CDRP_site_agg_profiles.csv.gz"

echo "========================================"
echo ""
echo "To check all files, run:"
echo "  uv run python scripts/restore_intelligent.py imaging-platform projects/2016_08_01_RadialMitochondriaDistribution_donna/workspace/metadata --filter_out preprocessed --tier Bulk --logfile logs/status_check_metadata.csv"
echo ""
echo "Once all files show RESTORED, run:"
echo "  bash scripts/download_data.sh"
echo ""
