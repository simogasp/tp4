# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added
- warning flags for gcc, clang and msvc
- added IPO support

### Changed
- modernize cmake
- switched to c++17

### Fixed
- all warnings on linux and windows

### Removed
- support for makefile building

## [2021.0.0]
### Added
- clang format

### Fixed
- fixed tests
- fixed generate student to copy only necessary files
- normal were shown wrong

### Changed
- better management of the input file with safeguard and messages

## [2020.1.0] 2020-05-14
### Added
- support for windows
- embed freeglut for windows
- ci on appveyor

### Changed
- flat structure of the code
- modernize the code, rely only on c++11
- raised cmake version


## [0.6.0] - 2017-10-05
First stable version with subdivision
### Added
- Subdivision


## [0.5.0] - 2014-10-19
First version of the visualizer

