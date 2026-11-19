"""Run the whole pipeline on SYNTHETIC data and stamp every output as a demo.

    python demo/run_demo.py      (about 4 minutes)

The facility data used in the real project are not public. This script invents a data set in the
same format (demo/make_synthetic_data.py), runs the same code on it in demo/workspace/, and marks
every figure and report as synthetic. The results here say nothing about any real system.
"""
import runpy
import shutil
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
WORK = HERE / "workspace"
MSG = "SYNTHETIC DEMO DATA - not real facility data"

if WORK.exists():
    shutil.rmtree(WORK)
shutil.copytree(HERE.parent / "src", WORK / "src")
sys.path.insert(0, str(HERE))
import make_synthetic_data  # noqa: E402

make_synthetic_data.main(WORK / "data" / "raw_anonymised")
sys.path.insert(0, str(WORK / "src"))
for script in ["audit.py", "conductivity_model.py", "ph_model.py", "calibrate_validate.py",
               "transfer_test.py", "ml_compare.py"]:
    print(f"\n=== {script} ===", flush=True)
    runpy.run_path(str(WORK / "src" / script), run_name="__main__")

# stamp the outputs
from matplotlib import font_manager  # noqa: E402
from PIL import Image, ImageDraw, ImageFont  # noqa: E402

FONT = ImageFont.truetype(font_manager.findfont("DejaVu Sans:bold"), 30)
BAR = 60


def stamp(png: Path) -> None:
    im = Image.open(png).convert("RGB")
    out = Image.new("RGB", (im.width, im.height + BAR), (170, 20, 20))
    ImageDraw.Draw(out).text((20, 12), MSG, fill=(255, 255, 255), font=FONT)
    out.paste(im, (0, BAR))
    out.save(png)


for png in (WORK / "figures").glob("*.png"):
    stamp(png)
for md in (WORK / "reports").glob("*.md"):
    md.write_text(f"> **{MSG}.** Results below describe invented data.\n\n" + md.read_text(encoding="utf-8"),
                  encoding="utf-8")
print("\nDone. See demo/workspace/figures and demo/workspace/reports.")
