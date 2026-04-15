from __future__ import annotations

import base64
import shutil
from pathlib import Path

import nbformat
from nbclient import NotebookClient


ROOT = Path(__file__).resolve().parents[1]
NOTEBOOK_PATH = ROOT / "README.ipynb"
README_PATH = ROOT / "README.md"
FIGURE_DIR = ROOT / "README_files" / "figure-markdown"


def ensure_clean_figure_dir() -> None:
    if FIGURE_DIR.parent.exists():
        shutil.rmtree(FIGURE_DIR.parent)
    FIGURE_DIR.mkdir(parents=True, exist_ok=True)


def save_png(image_b64: str, cell_index: int, output_index: int) -> str:
    filename = f"cell-{cell_index:02d}-output-{output_index:02d}.png"
    path = FIGURE_DIR / filename
    path.write_bytes(base64.b64decode(image_b64))
    return f"README_files/figure-markdown/{filename}"


def render_markdown(nb: nbformat.NotebookNode) -> str:
    chunks: list[str] = []
    for cell_index, cell in enumerate(nb.cells):
        if cell.cell_type == "markdown":
            chunks.append(cell.source.rstrip())
            continue

        if cell.cell_type != "code":
            continue

        hide_input = bool(cell.metadata.get("render", {}).get("hide_input"))
        if not hide_input and cell.source.strip():
            chunks.append(f"```python\n{cell.source.rstrip()}\n```")

        for output_index, output in enumerate(cell.get("outputs", [])):
            output_type = output.get("output_type")
            if output_type == "stream":
                text = output.get("text", "").rstrip()
                if text:
                    chunks.append(f"```text\n{text}\n```")
                continue

            data = output.get("data", {})
            if "image/png" in data:
                relpath = save_png(data["image/png"], cell_index, output_index)
                chunks.append(f"![]({relpath})")
                continue

            text_plain = data.get("text/plain")
            if text_plain:
                text = text_plain if isinstance(text_plain, str) else "".join(text_plain)
                chunks.append(f"```text\n{text.rstrip()}\n```")

        chunks.append("")
    return "\n\n".join(chunk for chunk in chunks if chunk is not None).rstrip() + "\n"


def main() -> None:
    nb = nbformat.read(NOTEBOOK_PATH, as_version=4)
    client = NotebookClient(nb, timeout=1200, kernel_name="python3", resources={"metadata": {"path": str(ROOT)}})
    client.execute()
    ensure_clean_figure_dir()
    README_PATH.write_text(render_markdown(nb), encoding="utf-8")


if __name__ == "__main__":
    main()
