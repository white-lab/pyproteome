# An example file showing batch searching of raw files using
# ProteomeDiscoverer 1.4


$PDD="C:\Program Files\Thermo\Discoverer 1.4\System\Release\DiscovererDaemon.exe"
$BASEDIR="Z:\morshed"
$WORKFLOWDIR="$BASEDIR\Workflows"


# Example ProteomeDiscoverer 1.4 Workflows are available online at:
# https://github.com/white-lab/pyproteome-data/tree/master/proteome_discoverer
$QE_6m_WORKFLOW="$WORKFLOWDIR\NFM_QE_dynTMT6plex_mouse_20mmu_pSTY.xml"
$QE_10m_WORKFLOW="$WORKFLOWDIR\NFM_QE_dynTMT10plex_mouse_20mmu_pSTY.xml"

$QE_6h_WORKFLOW="$WORKFLOWDIR\NFM_QE_dynTMT6plex_human_20mmu_pSTY.xml"
$QE_10h_WORKFLOW="$WORKFLOWDIR\NFM_QE_dynTMT10plex_human_20mmu_pSTY.xml"

$LTQ_6m_WORKFLOW="$WORKFLOWDIR\NFM_LTQ_Orbitrap_CIDHCD_mouse_dynTMT6plex_20ppm.xml"
$LTQ_10m_WORKFLOW="$WORKFLOWDIR\NFM_LTQ_Orbitrap_CIDHCD_mouse_dynTMT10plex_20ppm.xml"

$LTQ_6h_WORKFLOW="$WORKFLOWDIR\NFM_LTQ_Orbitrap_CIDHCD_human_dynTMT6plex_20ppm.xml"
$LTQ_610h_WORKFLOW="$WORKFLOWDIR\NFM_LTQ_Orbitrap_CIDHCD_human_dynTMT10plex_20ppm.xml"


# Run ProteomeDiscoverer Workflow on a set of Files, titled as Name.
# Save to DataDir\Name.msf
function run_pd($Name, $Files, $DataDir, $Workflow) {
    $msf_path = "$DataDir\$Name.msf"

    if ((Test-Path $msf_path) -and ((Get-Item $FileExists1).length -gt 1mb)) {
        return
    }

    $args = @()

    ForEach ($file in $Files) {
        if (!(Test-Path $file)) {
            echo "!!!!!!!!!!!"
            echo "! Warning ! $file is missing from your raw data collection!"
            echo "!!!!!!!!!!!"
        }
        $args += @("-a", $Name, $file)
    }

    $args += @(
        "-e", "$Name", $Files.Length, "$Workflow",
        "-r", "$msf_path"
    )

    & $PDD $args
}


$DATA_DIR="$BASEDIR\CK-p25"
$QE_DIR="$DATA_DIR\QE"
$LTQ_DIR="$DATA_DIR\LTQ"

$OUT_DIR=$DATA_DIR


run_pd "CK-H1-pY" @(
    "$QE_DIR\2015-09-11-CKH1-pY-imac14-elute-pre35-colAaron250.raw"
) $OUT_DIR $QE_6m_WORKFLOW
run_pd "CK-H1-MPM2" @(
    "$QE_DIR\2015-11-23-CKH1-MPM2-NTA-elute-pre43-col35.raw"
) $OUT_DIR $QE_6m_WORKFLOW
run_pd "CK-H1-pST" @(
    "$QE_DIR\2018-02-02-CKH1-MPM2-sup-pST-NoCond-pre140-col136.raw"
) $OUT_DIR $QE_6m_WORKFLOW
run_pd "CK-H1-Global" @(
    "$LTQ_DIR\2015-09-18-CKH1-pY-2-sup-10-preRaven-colAaron250.raw"
) $OUT_DIR $QE_6m_WORKFLOW


run_pd "CK-X11-pY" @(
    "$QE_DIR\2016-06-02-CKX11-pY-NTA-elute-pre60-col68.raw"
) $OUT_DIR $QE_6m_WORKFLOW
run_pd "CK-X11-MPM2" @(
    "$QE_DIR\2016-06-02-CKX11-MPM2-NTA-elute-pre63-col68.raw"
) $OUT_DIR $QE_6m_WORKFLOW
run_pd "CK-X11-pST" @(
    "$QE_DIR\2016-07-17-CKX11-pST-1-NTA-elute.raw",
    "$QE_DIR\2016-07-17-CKX11-pST-2-NTA-elute.raw",
    "$QE_DIR\2016-07-17-CKX11-pST-3-NTA-elute.raw",
    "$QE_DIR\2016-07-17-CKX11-pST-4-NTA-elute.raw",
    "$QE_DIR\2016-07-17-CKX11-pST-5-NTA-elute.raw",
    "$QE_DIR\2016-07-17-CKX11-pST-6-NTA-elute.raw",
    "$QE_DIR\2016-07-17-CKX11-pST-7-NTA-elute.raw",
    "$QE_DIR\2016-07-17-CKX11-pST-8-NTA-elute.raw",
    "$QE_DIR\2016-07-17-CKX11-pST-10-NTA-elute.raw",
    "$QE_DIR\2016-07-17-CKX11-pST-11-NTA-elute.raw",
    "$QE_DIR\2016-07-17-CKX11-pST-15-NTA-elute.raw",
    "$QE_DIR\2016-07-17-CKX11-pST-16-NTA-elute.raw",
    "$QE_DIR\2016-07-17-CKX11-pST-17-NTA-elute.raw",
    "$QE_DIR\2016-07-17-CKX11-pST-18-NTA-elute.raw",
    "$QE_DIR\2016-07-17-CKX11-pST-19-NTA-elute.raw",
    "$QE_DIR\2016-07-17-CKX11-pST-20-NTA-elute.raw"
) $OUT_DIR $QE_6m_WORKFLOW
run_pd "CK-X11-Global" @(
    "$LTQ_DIR\2016-06-02-CKX11-pY-sup10-pre63-col68.raw"
) $OUT_DIR $LTQ_6m_WORKFLOW
