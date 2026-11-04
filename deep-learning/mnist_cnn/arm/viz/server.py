"""Dashboard server, Python standard library only. Reads the files that mnistarm writes."""
import json
import re
import subprocess
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path

ARM = Path(__file__).resolve().parent.parent
RUNS = ARM / "runs"
BINARY = ARM / "build" / "mnistarm"
PORT = 8577
bench_cache = {}


def read_json(path):
    try:
        return json.loads(path.read_text())
    except (OSError, ValueError):
        return None


def run_data(name):
    run = RUNS / name
    metrics = []
    if (run / "metrics.jsonl").exists():
        for line in (run / "metrics.jsonl").read_text().splitlines():
            try:
                metrics.append(json.loads(line))
            except ValueError:
                pass
    config = {}
    if (run / "config.txt").exists():
        config = dict(line.split() for line in (run / "config.txt").read_text().splitlines() if line)
    return {"name": name, "metrics": metrics, "config": config,
            "snapshot": read_json(run / "snapshot.json"), "test": read_json(run / "test.json")}


class Handler(BaseHTTPRequestHandler):
    def log_message(self, *args):
        pass

    def send(self, body, kind="application/json", code=200):
        data = body if isinstance(body, bytes) else json.dumps(body).encode()
        self.send_response(code)
        self.send_header("Content-Type", kind)
        self.send_header("Content-Length", str(len(data)))
        self.send_header("Cache-Control", "no-store")
        self.end_headers()
        self.wfile.write(data)

    def do_GET(self):
        if self.path == "/":
            return self.send((ARM / "viz" / "index.html").read_bytes(), "text/html; charset=utf-8")
        if self.path == "/api/runs":
            names = sorted(p.parent.name for p in RUNS.glob("*/metrics.jsonl"))
            names.sort(key=lambda n: (RUNS / n / "metrics.jsonl").stat().st_mtime, reverse=True)
            return self.send(names)
        if self.path == "/api/bench":
            if not bench_cache:
                out = subprocess.run([BINARY, "bench"], cwd=ARM, capture_output=True, text=True)
                bench_cache.update(json.loads(out.stdout))
            return self.send(bench_cache)
        match = re.fullmatch(r"/api/run/([\w.-]+)", self.path)
        if match and (RUNS / match[1]).is_dir():
            return self.send(run_data(match[1]))
        self.send({"error": "not found"}, code=404)

    def do_POST(self):
        match = re.fullmatch(r"/api/predict/([\w.-]+)", self.path)
        if not match or not (RUNS / match[1] / "model.bin").exists():
            return self.send({"error": "no such run"}, code=404)
        pixels = json.loads(self.rfile.read(int(self.headers["Content-Length"])))
        if len(pixels) != 784:
            return self.send({"error": "need 784 values"}, code=400)
        text = " ".join(f"{float(v):.4f}" for v in pixels)
        out = subprocess.run([BINARY, "predict", "--run", match[1]], cwd=ARM, input=text,
                             capture_output=True, text=True)
        self.send(out.stdout.encode() if out.returncode == 0 else b'{"error":"predict failed"}')


if __name__ == "__main__":
    print(f"http://localhost:{PORT}")
    ThreadingHTTPServer(("127.0.0.1", PORT), Handler).serve_forever()
