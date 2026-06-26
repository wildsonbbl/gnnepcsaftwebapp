param(
	[switch]$SkipUpload
)

$versionFile = Join-Path $PSScriptRoot 'webapp/_version.py'
$versionNumber = Select-String -Path $versionFile -Pattern '"([0-9]+\.[0-9]+\.[0-9]+)"' |
	Select-Object -First 1 |
	ForEach-Object { $_.Matches[0].Groups[1].Value }

if (-not $versionNumber) {
	throw "Could not find version in $versionFile"
}

$version = "v$versionNumber"
$platform='windows'
$installerName = "gnnpcsaftwebapp-$version-$platform.msi"

## create package
uv pip install -r requirements.txt
uv run python manage.py collectstatic --no-input
uv run python manage.py migrate --no-input
uv run pyinstaller --distpath ./app_pkg/dist --workpath ./app_pkg/build --noconfirm --clean ./gnnpcsaftwebapp.spec

$distDir = Join-Path $PSScriptRoot 'app_pkg/dist/gnnpcsaftwebapp'
if (-not (Test-Path $distDir)) {
	throw "Could not find app_pkg/dist directory at $distDir"
}

## create installer
$installerOutputDir = Join-Path $PSScriptRoot 'installer'
$productWxs = Join-Path $PSScriptRoot 'gnnpcsaftwebapp-product.wxs'

Remove-Item -Path $installerOutputDir -Recurse -Force -ErrorAction SilentlyContinue
New-Item -ItemType Directory -Path $installerOutputDir -Force | Out-Null

$installerArtifact = Join-Path $installerOutputDir $installerName
$wixCommand = Get-Command wix.exe -ErrorAction SilentlyContinue | Select-Object -First 1

if ($wixCommand) {
	$wixExe = if ($wixCommand.Source) { $wixCommand.Source } else { $wixCommand.Path }
	& $wixExe build --acceptEula wix7 $productWxs -arch x64 -d "ProductVersion=$versionNumber" -d "SourceDir=$distDir" -d "ProjectDir=$PSScriptRoot" -o $installerArtifact
	if ($LASTEXITCODE -ne 0) {
		throw 'WiX build failed while generating MSI with wix.exe'
	}
}

if (-not (Test-Path $installerArtifact)) {
	throw "Could not find generated installer at $installerArtifact"
}

## add artifact to release
if (-not $SkipUpload) {
	gh release upload $version $installerArtifact
}