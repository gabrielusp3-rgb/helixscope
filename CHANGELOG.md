# Changelog

The application version recorded in results is `0.24.3-19`.
Earlier development notes in this repository are working records, not a
reconstructed release history.

## Unreleased

### Added

- Linux container for the Streamlit workstation, with FastTree 2.1.11,
  IQ-TREE 3.1.3, NCBI BLAST+ 2.17.0+, US-align 20260328 and ViennaRNA 2.7.2.
- Runtime `PORT` selection. Absent `PORT` binds 8501. Invalid values stop
  startup. The health check uses the same port.
- Session analysis history, with separate entries per browser session.

### Changed

- The public README keeps the repository banner and no longer embeds the
  interface screenshots or the captions that described them.
- Python dependencies in `requirements.txt` are pinned to the versions measured
  in the tested container: Streamlit 1.58.0, Biopython 1.87, Plotly 6.9.0,
  pandas 2.3.3, NumPy 2.5.3, pytest 9.1.1 and Hypothesis 6.168.3.
- The container base is Debian Trixie with Python 3.14.3, after Bookworm
  security updates left two critical Perl issues without a fix.

### Fixed

- Windows tool paths are recognized on Linux.
- The first navigation away from Overview can return to Overview.
- The container user can use a read-only root filesystem when home and `/tmp`
  are mounted as writable temporary directories.

### Security

- The image runs as `helix` (uid 10001), with a read-only root in Compose,
  a 64-process limit and a 1 GiB memory limit.
- Docker Scout on the tested image reported no critical OS vulnerabilities.
  Three high findings remained in `gcc-14` runtime libraries and `zlib`, with
  no Debian fix at scan time.
- Optional engine absence stays an explicit not-installed state.

## [0.24.3-19]

Declared provenance version. This file does not reconstruct the changes that
led to that number.
