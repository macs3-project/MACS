MACS3.IO.Parser module
======================

BAM/BAMPE and FRAG parsers share an internal compiled byte stream in
``MACS3.IO.BAMStream``. Public parser imports and methods remain in
``MACS3.IO.Parser``. The stream owns decompression, native workers and
compressed/decompressed buffers; its internal Cython declarations describe
the buffer lifetime and serial-fallback contracts. Parser retains format
decoding, FRAG interning, fragment batches and track assembly. Extracting
those FRAG components is deferred because they share native callback layouts
and track emission logic.

.. automodule:: MACS3.IO.Parser
   :members:
   :undoc-members:
   :show-inheritance:
   :exclude-members: d, n, StrandFormatError, atoi, bool,bam_fw_binary_parse, bampe_pe_binary_parse, tlen_parse_line
