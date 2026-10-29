"""Validation shared by commands that require alignment tracks."""


def require_nonempty_track(track, input_format, role="input"):
    """Raise a descriptive error when loading or filtering kept no alignments."""
    if track.total == 0:
        message = f"No usable alignments remain in {role} after loading or filtering."
        if input_format in ("BAM", "BAMPE"):
            message += " Paired BAM reads require proper pairs (flag 0x2)."
        raise ValueError(message)
