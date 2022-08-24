for %%f in (./inputs/*) do (
  call a.exe < ./inputs/%%f > ./inputs2/%%f
)