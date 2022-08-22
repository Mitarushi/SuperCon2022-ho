for %%f in (./inputs/*) do (
  call a.exe < ./inputs/%%f > ./outputs/%%f
)