from PyInstaller.utils.hooks import collect_submodules

# openseespy resolves solver imports dynamically at runtime. 
# We explicitly collect all submodules of both openseespy and openseespywin.
hiddenimports = collect_submodules('openseespy')
hiddenimports += collect_submodules('openseespywin')
