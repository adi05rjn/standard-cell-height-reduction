# Points a translated .vproc at the project's material file (the JSON Translator leaves this blank).
# Usage: python link_materials.py <folder> <name>
import io, os, re, sys
folder, name = sys.argv[1], sys.argv[2]
p = os.path.join(folder, name + ".vproc")
s = io.open(p, encoding="utf-8").read()
absp = os.path.abspath(os.path.join(folder, name + ".vmpd")).replace("\\", "/")
s2, n = re.subn(r'<materialDatabase RelativeFileName="[^"]*" ___text="[^"]*"',
                '<materialDatabase RelativeFileName="./%s.vmpd" ___text="%s"' % (name, absp), s, count=1)
io.open(p, "w", encoding="utf-8", newline="\n").write(s2)
print("material link set:", n == 1, "steps:", s2.count("<stepName"))
