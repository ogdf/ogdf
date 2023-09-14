# Runs concurrent test on windows.
#
# Author: Tilo Wiedera

param ([switch]$debug, [switch]$gurobi)

$modeFlag = "release"
$dir = "Release"
If ( $debug ) {
	$modeFlag = "debug"
	$dir = "Debug"
}
$env:path = "$env:path;$pwd\$dir"

$ilpsolver = "default_s"
If ( $gurobi ) {
    $ilpsolver = "gurobi"
}

& .\util\build.bat $modeFlag $ilpsolver
If ( $LASTEXITCODE -ne 0 ) {
	Exit 1
}

$id = Get-Random
$pattern = "$id *"

Get-ChildItem "test\bin\$dir" -Filter *.exe | Foreach-Object {
	$ScriptBlock = {
		param($filename)
		$testlog = & "test\bin\$using:dir\$filename" --break-on-failure --reporter=crash 2>&1
		If (-Not $?) {
			$errors = Select-String -InputObject $testlog -Pattern "^#" -Context 1,5
			If (-Not $errors) {
				$errors = Write-Output -InputObject $testlog | Select-Object -Last 10
			}
			$errors = $errors -join "`r`n"
			Write-Error -ErrorAction Stop -Message "Test Error: $filename`r`n$errors`r`n"
		}
	}

	$name = $_.BaseName
	Start-Job -Name "$id $name" -Init ([ScriptBlock]::Create("Set-Location $pwd")) -Scrip $ScriptBlock -ArgumentList $_ | select State, Name
}

Get-ChildItem "doc\examples" -Recurse -Filter ex-*.exe -Exclude ex-multilevelmixer.exe | Foreach-Object {
	$ScriptBlock = {
		param($filename)
		& "$filename"
		If (-Not $?) {
			Write-Error -ErrorAction Stop -Message "Example Error: $filename`r`n"
		}
	}

	$name = $_.BaseName
	$dir = $_.DirectoryName
	Start-Job -Name "$id $name" -Init ([ScriptBlock]::Create("Set-Location $dir\..")) -Scrip $ScriptBlock -ArgumentList $_ | select State, Name
}

$errorCode = 0

function Check-Failure {
	$failedJobs = Get-Job -Name $pattern | Where-Object { $_.State -eq "Failed" }
	If ($failedJobs) {
		$failedJobs | Receive-Job
		$script:errorCode = 1
	}
	return $failedJobs
}

While (Get-Job -Name $pattern | Where-Object { $_.State -eq "Running" })
{
	If (Check-Failure) {
		break
	}
	Start-Sleep 1
}
Check-Failure | Out-Null

Get-Job -Name $pattern | select State, Name
Get-Job -Name $pattern | Remove-Job -Force

Exit $errorCode
