; Windows installer for the tdtK Heart Analyzer (Inno Setup 6).
; Built by the GitHub workflow after PyInstaller:
;   iscc /DAppVersion=1.0.0a1 packaging\installer.iss
; Input: packaging\dist\tdtK-Heart-Analyzer\   Output: packaging\installer\tdtK-Heart-Analyzer-Setup-<version>.exe
;
; No administrator rights needed: by default it installs for the current user only;
; the first page offers "install for all users" (needs admin) instead.

#ifndef AppVersion
  #define AppVersion "0.0.0"
#endif
#define AppName "tdtK Heart Analyzer"
#define AppExe "tdtK-Heart-Analyzer.exe"

[Setup]
AppId={{8C4018FB-3FFD-4195-80EA-DE9BE2491B1D}
AppName={#AppName}
AppVersion={#AppVersion}
AppVerName={#AppName} {#AppVersion}
AppPublisher=FlyHearts
AppPublisherURL=https://github.com/gvogler/FlyHearts-tdtK-Rscripts
AppSupportURL=https://github.com/gvogler/FlyHearts-tdtK-Rscripts
DefaultDirName={autopf}\{#AppName}
DefaultGroupName={#AppName}
DisableProgramGroupPage=yes
PrivilegesRequired=lowest
PrivilegesRequiredOverridesAllowed=dialog
ArchitecturesAllowed=x64compatible
ArchitecturesInstallIn64BitMode=x64compatible
LicenseFile=..\..\LICENSE
SetupIconFile=tdtk-analyzer.ico
UninstallDisplayIcon={app}\{#AppExe}
OutputDir=installer
OutputBaseFilename=tdtK-Heart-Analyzer-Setup-{#AppVersion}
Compression=lzma2/max
SolidCompression=yes
WizardStyle=modern

[Languages]
Name: "english"; MessagesFile: "compiler:Default.isl"

[Tasks]
Name: "desktopicon"; Description: "{cm:CreateDesktopIcon}"; GroupDescription: "{cm:AdditionalIcons}"

[Files]
Source: "dist\tdtK-Heart-Analyzer\*"; DestDir: "{app}"; Flags: ignoreversion recursesubdirs createallsubdirs

[Icons]
Name: "{autoprograms}\{#AppName}"; Filename: "{app}\{#AppExe}"
Name: "{autodesktop}\{#AppName}"; Filename: "{app}\{#AppExe}"; Tasks: desktopicon

[Run]
Filename: "{app}\{#AppExe}"; Description: "{cm:LaunchProgram,{#AppName}}"; Flags: nowait postinstall skipifsilent
