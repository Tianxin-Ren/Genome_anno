import sys
import base64
import re
import os
import yaml

def calc_threads(requested, config, workflow=None):
    """Global helper to cap threads based on profile/config max_threads."""
    # Resolve max_threads from config or profile
    mt = config.get("max_threads")
    if mt is None and workflow is not None:
        _rs = getattr(workflow, "resource_settings", None)
        if _rs:
            default_res = getattr(_rs, "default_resources", None)
            if default_res is not None:
                if hasattr(default_res, "get"):
                    mt = default_res.get("max_threads")
                else:
                    mt = getattr(default_res, "max_threads", None)
        if mt is None:
            default_res = getattr(workflow, "default_resources", None)
            if default_res is not None:
                if hasattr(default_res, "get"):
                    mt = default_res.get("max_threads")
                else:
                    mt = getattr(default_res, "max_threads", None)
    
    if mt is not None:
        try:
            mt_val = int(getattr(mt, "value", mt))
            return min(int(requested), mt_val)
        except:
            pass
    return int(requested)

def calc_mem_mb(threads, workflow=None):
    """Global helper to calculate memory based on profile/config default_resources mem_mb."""
    multiplier = 8000 # fallback
    
    # 1. Check sys.argv for --default-resources (Used by snakemake remote executor)
    in_default_res = False
    for arg in sys.argv:
        if arg == "--default-resources":
            in_default_res = True
            continue
        if arg.startswith("-"):
            in_default_res = False
            
        if in_default_res:
            val = arg
            if arg.startswith("base64//"):
                try:
                    val = base64.b64decode(arg[8:]).decode("utf-8")
                except:
                    pass
            if val.startswith("mem_mb="):
                m = re.search(r'(\d+)', val)
                if m:
                    return threads * int(m.group(1))
                    
    # 2. Check sys.argv for --profile or --workflow-profile and read its config.yaml
    for i, arg in enumerate(sys.argv):
        if arg in ("--profile", "--workflow-profile") and i + 1 < len(sys.argv):
            prof_path = sys.argv[i+1]
            conf_file = os.path.join(prof_path, "config.yaml")
            if os.path.exists(conf_file):
                try:
                    with open(conf_file) as f:
                        p_conf = yaml.safe_load(f)
                        mem_str = p_conf.get("default-resources", {}).get("mem_mb", "")
                        m = re.search(r'(\d+)', str(mem_str))
                        if m:
                            return threads * int(m.group(1))
                except:
                    pass

    # 3. Check workflow default_resources object (fallback)
    if workflow is not None:
        _rs = getattr(workflow, "resource_settings", None)
        mem_str = None
        if _rs:
            default_res = getattr(_rs, "default_resources", None)
            if default_res is not None:
                if hasattr(default_res, "get"):
                    mem_str = default_res.get("mem_mb")
                else:
                    mem_str = getattr(default_res, "mem_mb", None)
        if mem_str is None:
            default_res = getattr(workflow, "default_resources", None)
            if default_res is not None:
                if hasattr(default_res, "get"):
                    mem_str = default_res.get("mem_mb")
                else:
                    mem_str = getattr(default_res, "mem_mb", None)
            
        if mem_str is not None:
            try:
                val = str(getattr(mem_str, "value", mem_str))
                m = re.search(r'(\d+)', val)
                if m:
                    return threads * int(m.group(1))
            except:
                pass
            
    return threads * multiplier
