#!/bin/bash
# Start the Open Meta Maps tool and open it in the browser.
# The server stops by itself when the page is closed (and no job is running); closing this window also stops it.
cd "$(dirname "$0")" || exit 1
PORT="${1:-8765}"
export PATH="$HOME/.local/node/bin:$PATH" PYTHONHASHSEED=0

# Python with the drawing dependencies (PIL, jsonschema, scipy, networkx, numpy): ~/.local/opm-venv (outside the synced folder), then a project .venv, the old test venv, then system python
for PY in "$HOME/.local/opm-venv/bin/python3" "$PWD/.venv/bin/python3" "/private/tmp/claude-501/venv2/bin/python3" "$(command -v python3)"; do
  if [ -x "$PY" ] && "$PY" -c "import PIL, jsonschema, scipy, networkx, numpy" 2>/dev/null; then break; fi
  PY=""
done
if [ -z "$PY" ]; then
  echo "没有找到装有 PIL、jsonschema、scipy、networkx、numpy 的 Python。"
  echo "Python with PIL, jsonschema, scipy, networkx and numpy was not found."
  read -r -p "按回车关闭 / press Return to close " _; exit 1
fi

echo "Python: $PY"
"$PY" layout_atlas/cli.py tool --port "$PORT" &
SERVER=$!
trap 'kill $SERVER 2>/dev/null' EXIT INT TERM HUP

# wait until the server answers, then open the page (the server may pick another port if $PORT belongs to a different program)
for i in $(seq 1 40); do
  sleep 0.5
  kill -0 $SERVER 2>/dev/null || { echo "服务没有启动 / the server did not start"; read -r -p "按回车关闭 / press Return to close " _; exit 1; }
  for P in $(seq "$PORT" $((PORT + 29))); do
    if curl -s -m 1 "http://127.0.0.1:$P/api/health" | grep -q open-pathway-maps; then
      open "http://127.0.0.1:$P/tool/index.html"; break 2
    fi
  done
done
wait $SERVER
echo "工具已停止 / the tool has stopped."
