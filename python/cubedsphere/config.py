"""配置解析（YAML）。

生产环境优先使用 PyYAML；未安装时回退到内置的极简 YAML 子集解析器
（支持嵌套映射、标量、注释），覆盖本项目的 default.yaml 结构。
"""

from __future__ import annotations


def load_yaml(text: str) -> dict:
    """解析 YAML 文本为嵌套 dict。

    优先 PyYAML；回退到内置子集解析器。
    """
    try:
        import yaml  # type: ignore
    except ImportError:
        return _parse_simple_yaml(text)
    return yaml.safe_load(text)


def load_yaml_file(path: str) -> dict:
    with open(path, "r", encoding="utf-8") as f:
        return load_yaml(f.read())


def _parse_simple_yaml(text: str) -> dict:
    """极简 YAML 子集解析器（嵌套映射 + 标量 + 注释）。

    支持：缩进嵌套、``key: value`` 标量、整数/浮点/布尔/字符串、
    行内 ``#`` 注释。不支持列表、块引用、多行字符串等高级特性。
    """
    result: dict = {}
    # 栈保存 (缩进层级, dict 引用)
    stack: list[tuple[int, dict]] = [(-1, result)]

    for raw_line in text.splitlines():
        line = raw_line.split("#")[0].rstrip()
        if not line.strip():
            continue
        indent = len(line) - len(line.lstrip(" "))
        content = line.strip()
        if ":" not in content:
            continue
        key, _, value = content.partition(":")
        key = key.strip()
        value = value.strip()

        # 弹栈到合适缩进
        while stack and indent <= stack[-1][0]:
            stack.pop()
        parent = stack[-1][1]

        if value == "":
            # 嵌套映射开始
            child: dict = {}
            parent[key] = child
            stack.append((indent, child))
        else:
            parent[key] = _parse_scalar(value)

    return result


def _parse_scalar(s: str):
    """解析 YAML 标量为 Python 原生类型。"""
    low = s.lower()
    if low == "true":
        return True
    if low == "false":
        return False
    if low in ("null", "~"):
        return None
    # 去掉引号
    if len(s) >= 2 and s[0] == s[-1] and s[0] in "\"'":
        return s[1:-1]
    try:
        return int(s)
    except ValueError:
        pass
    try:
        return float(s)
    except ValueError:
        pass
    return s


def get_path(config: dict, dotted_key: str, default=None):
    """按点号路径读取配置值（grid.ncells → config["grid"]["ncells"]）。"""
    cur = config
    for part in dotted_key.split("."):
        if not isinstance(cur, dict) or part not in cur:
            return default
        cur = cur[part]
    return cur
