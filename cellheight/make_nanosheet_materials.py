# Creates common/cellheight_nanosheet.vmpd: the 16 nm example materials plus "SiGe_Sacrificial",
# a renamed copy of SiGe, so the sacrificial layers can be removed without touching the P source/drain.
# Usage (inside the SEMulator3D shell): python make_nanosheet_materials.py
_SEMulator3D_API_VERSION = '13.1'
import sys, os
from voxelModeler import *
from CModeler import mpiInitialize, mpiFinalize
mpiInitialize(len(sys.argv), sys.argv)
here = os.path.dirname(os.path.abspath(__file__))
matDoc = CMaterialData.MaterialDocument()
matDoc.openWithPasswords(os.path.join(here, "common", "16nmFINFET.vmpd"), "", "", False)
mats = list(matDoc.getAllMaterials())
sige = [m for m in mats if m.getName() == "SiGe"][0] if hasattr(mats[0], "getName") else None
print("MATOBJ", [n for n in dir(mats[0]) if not n.startswith("_")])
print("SIGE", sige.getName(), sige.getFullUntypedName() if hasattr(sige, "getFullUntypedName") else "")
root = matDoc.getMaterialRoot()
print("ROOTCHILDREN", root.getChildCount())
grp = None
for i in range(root.getChildCount()):
    try:
        c = root.getChild(i)
        print("  child", i, c.getName(), c.getType())
        if c.getName() == "Silicon": grp = c
    except Exception as e:
        print("  getChild failed", e); break
try:
    new = sige.clone()
    new.setName("SiGe_Sacrificial")
    try: new.newUUID()
    except Exception as e: print("newUUID", e)
    (grp or root).addChild(new)
    out = os.path.join(here, "common", "cellheight_nanosheet.vmpd")
    print("SAVE", matDoc.save(out) if True else None, os.path.exists(out))
    m2 = CMaterialData.MaterialDocument(); m2.openWithPasswords(out, "", "", False)
    print("REOPEN", [m.getFullUntypedName() for m in m2.getAllMaterials() if "SiGe" in m.getName()])
except Exception as e:
    import traceback; traceback.print_exc()
mpiFinalize()
