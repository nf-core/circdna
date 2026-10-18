# Hotfix for quay.io/bioinfortools/fled:1.7.0 (built 2026-09-20).
#
# The image was built from a context containing a stale setuptools `build/`
# directory. `pip install .` reused that artifact instead of the current
# source, so the installed copy of FLED/longreadscircDNAdetection.py lost the
# `import os` line while still calling os.getpid()/os.getcwd(). FLED therefore
# died at startup with:
#
#     NameError: name 'os' is not defined
#
# The stale host artifact has since been synced with the source, so a full
# rebuild from /data1/users/siyangming/opt/FLED is also correct. This layer
# makes the already-built tag usable without a multi-minute rebuild.
FROM quay.io/bioinfortools/fled:1.7.0

RUN set -eu; \
    f=/usr/local/lib/python3.8/site-packages/FLED/longreadscircDNAdetection.py; \
    sed -i '1a import os' "$f"; \
    grep -q '^import os' "$f"; \
    b=/opt/FLED/build/lib/FLED/longreadscircDNAdetection.py; \
    if [ -f "$b" ]; then sed -i '1a import os' "$b"; fi

# The bug only surfaced at runtime (os.getcwd() inside args_detect), so check
# that the module actually has `os` bound rather than merely that it imports.
RUN python -c "import FLED.longreadscircDNAdetection as m; assert hasattr(m, 'os'), 'os missing'; print('FLED os import OK:', m.os.getcwd())"
