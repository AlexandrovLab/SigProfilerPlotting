# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/), and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

## [1.4.1] - 2025-03-10

### Fixed
- The plotting function for structural variants has had syntax to reflect pandas 2.0.0 changes.

## [1.4.0] - 2025-02-11

### Changed
- Updated dependencies: Now requires **Pandas >= 2.0.0**, **NumPy >= 2.0.0**, and **Python >= 3.9**.
- Dropped support for **Python 3.8**

## [1.3.24] - 2024-08-01

### Added
- Added environmental variable `SIGPROFILERPLOTTING_VOLUME` to enhance configurability.

### Fixed
- Ensured correct file path construction with `os.path.join()`.
- Updated CLI to handle boolean arguments using `str2bool`, resolving an issue where boolean strings like `False` were incorrectly evaluated as `True`.
- Stopped unecessary output directory creation when the output format is PIL_Image.
