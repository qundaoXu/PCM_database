#!/usr/bin/env python3
"""Start or reuse the read-only localhost preview for docs/ in pcm_database."""
import argparse
from functools import partial
from http.server import SimpleHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path
from urllib.error import URLError
from urllib.request import urlopen
import webbrowser

DOCS = Path(__file__).resolve().parents[1] / "docs"


def is_catalog_preview(url):
    try:
        with urlopen(url + "index.html", timeout=2) as response:
            html = response.read().decode("utf-8")
        return 'id="tbl"' in html and 'src="catalog_logic.js"' in html
    except (URLError, OSError, UnicodeError):
        return False


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--port", type=int, default=8765)
    parser.add_argument("--no-browser", action="store_true", help="Serve or check without opening a browser")
    args = parser.parse_args()
    url = f"http://127.0.0.1:{args.port}/"
    if is_catalog_preview(url):
        print(f"数据库预览已在运行：{url}", flush=True)
        if not args.no_browser:
            webbrowser.open(url)
        return
    handler = partial(SimpleHTTPRequestHandler, directory=str(DOCS))
    try:
        server = ThreadingHTTPServer(("127.0.0.1", args.port), handler)
    except OSError as error:
        raise SystemExit(f"无法启动预览：端口 {args.port} 已占用或不可访问。{error}") from error
    with server:
        print(f"数据库预览已启动：{url}\n保持此窗口打开；按 Ctrl+C 关闭预览。", flush=True)
        if not args.no_browser:
            webbrowser.open(url)
        try:
            server.serve_forever()
        except KeyboardInterrupt:
            print("\n数据库预览已关闭。", flush=True)


if __name__ == "__main__":
    main()
