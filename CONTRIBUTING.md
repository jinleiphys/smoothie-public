# Contributing to SMOOTHIE

🎉 Thank you for your interest in contributing to SMOOTHIE! We welcome contributions from the community and are excited to work with you.

## 📋 Table of Contents

- [Getting Started](#getting-started)
- [Development Setup](#development-setup)
- [Contributing Guidelines](#contributing-guidelines)
- [Code Style](#code-style)
- [Testing](#testing)
- [Submitting Changes](#submitting-changes)
- [Issue Reporting](#issue-reporting)
- [Community Guidelines](#community-guidelines)

## 🚀 Getting Started

### Prerequisites

Before contributing, ensure you have:
- **Fortran 95/2003** compiler (gfortran 4.8+, ifort 14+)
- **Git** for version control
- **Make** build system
- **LAPACK/BLAS** libraries (recommended)
- Basic knowledge of quantum scattering theory (helpful but not required)

### Development Setup

1. **Fork the repository** on GitHub
2. **Clone your fork** locally:
   ```bash
   git clone https://github.com/your-username/smoothie-public.git
   cd smoothie-public
   ```

3. **Set up upstream remote**:
   ```bash
   git remote add upstream https://github.com/jinleiphys/smoothie-public.git
   ```

4. **Configure build environment**:
   ```bash
   cp make.inc.example make.inc
   vim make.inc  # Edit compiler settings
   ```

5. **Build and test**:
   ```bash
   cd smoothie/
   make
   cd test/
   ./smoothie < test.in
   ```

## 🤝 Contributing Guidelines

### Types of Contributions

We welcome several types of contributions:

#### 🐛 Bug Reports
- Found a bug? Please report it!
- Check existing issues first
- Include minimal reproducible example
- Provide system information

#### 🚀 Feature Requests
- Suggest new features or enhancements
- Explain the use case and benefits
- Discuss implementation approach

#### 📚 Documentation
- Improve existing documentation
- Add examples and tutorials
- Fix typos and clarify explanations

#### 🔧 Code Contributions
- Bug fixes
- Performance improvements
- New features
- Test additions

#### 📊 Input Examples
- Real-world calculation examples
- Validation cases
- Benchmark problems

### What We're Looking For

**High Priority:**
- 🔬 **Validation cases** and benchmark problems
- 📝 **Documentation improvements** and examples
- 🐛 **Bug fixes** and stability improvements
- ⚡ **Performance optimizations**

**Medium Priority:**
- 🆕 **New potential models** (with proper references)
- 🧪 **Additional test cases**
- 🔧 **Build system improvements**

**Please Discuss First:**
- 🏗️ **Major architectural changes**
- 📐 **New calculation methods**
- 🔄 **API changes**

## 📝 Code Style

### Fortran Coding Standards

**General Principles:**
- Follow **Fortran 95/2003** standards
- Use **meaningful variable names**
- Add **comments** for complex physics
- Maintain **consistent indentation** (2 spaces)

**Naming Conventions:**
```fortran
! Modules: lowercase with underscores
module quantum_scattering

! Variables: lowercase with underscores
real(dp) :: binding_energy
integer :: angular_momentum_max

! Constants: uppercase with underscores
real(dp), parameter :: PI = 3.14159265359_dp
```

**Documentation:**
```fortran
!> Brief description of subroutine
!> @param[in] energy - incident energy in MeV
!> @param[out] cross_section - calculated cross section
subroutine calculate_cross_section(energy, cross_section)
```

### File Organization

**Source Files:**
- `*.f90` - Main Fortran source files
- `*.f` - Legacy Fortran 77 files (minimize use)
- `*.inc` - Include files for common parameters

**Directory Structure:**
```
smoothie/
├── src/               # Main source code
├── modules/           # Fortran modules
├── test/             # Test cases and examples
├── doc/              # Documentation
└── examples/         # Example input files
```

## 🧪 Testing

### Test Requirements

**All contributions must:**
- ✅ **Compile successfully** on major platforms
- ✅ **Pass existing tests** without regression
- ✅ **Include new tests** for new features
- ✅ **Validate physics** results when applicable

### Running Tests

```bash
# Basic compilation test
make clean && make

# Run test suite
cd test/
./run_tests.sh

# Specific test case
./smoothie < d93nb.in
```

### Adding Tests

**For new features:**
1. Create input file in `test/`
2. Add expected output reference
3. Include test in `run_tests.sh`
4. Document the test case

**For bug fixes:**
1. Create minimal test case reproducing the bug
2. Verify fix resolves the issue
3. Add regression test

## 📤 Submitting Changes

### Pull Request Process

1. **Create a feature branch**:
   ```bash
   git checkout -b feature/your-feature-name
   ```

2. **Make your changes**:
   - Follow code style guidelines
   - Add tests for new functionality
   - Update documentation as needed

3. **Commit your changes**:
   ```bash
   git add .
   git commit -m "Add feature: brief description
   
   - Detailed description of changes
   - Why the change is needed
   - Any breaking changes"
   ```

4. **Push to your fork**:
   ```bash
   git push origin feature/your-feature-name
   ```

5. **Open a Pull Request**:
   - Use the PR template
   - Reference related issues
   - Describe changes and testing

### Pull Request Template

```markdown
## Description
Brief description of changes

## Type of Change
- [ ] Bug fix
- [ ] New feature
- [ ] Documentation update
- [ ] Performance improvement
- [ ] Other (please describe)

## Testing
- [ ] All tests pass
- [ ] New tests added
- [ ] Manual testing completed

## Checklist
- [ ] Code follows style guidelines
- [ ] Self-review completed
- [ ] Documentation updated
- [ ] No breaking changes (or clearly documented)
```

## 🐛 Issue Reporting

### Bug Reports

**Good bug reports include:**
- 🔍 **Clear description** of the problem
- 📋 **Steps to reproduce** the issue
- 🖥️ **Environment details** (OS, compiler, versions)
- 📁 **Input files** that demonstrate the problem
- 📊 **Expected vs actual behavior**

**Bug Report Template:**
```markdown
## Bug Description
Clear description of what the bug is.

## To Reproduce
Steps to reproduce the behavior:
1. Compile with '...'
2. Run with input '...'
3. See error

## Expected Behavior
What you expected to happen.

## Environment
- OS: [e.g., Ubuntu 20.04]
- Compiler: [e.g., gfortran 9.3]
- SMOOTHIE version: [e.g., commit hash]

## Additional Context
Any other context about the problem.
```

### Feature Requests

**Include:**
- 📝 **Use case description**
- 🎯 **Proposed solution**
- 🔬 **Physics background** (if applicable)
- 📚 **References** to literature

## 🌟 Community Guidelines

### Code of Conduct

We are committed to providing a welcoming and inclusive environment:

- 🤝 **Be respectful** and considerate
- 💬 **Communicate constructively**
- 🎓 **Help newcomers** learn
- 📚 **Share knowledge** freely
- 🔬 **Focus on physics** and code quality

### Getting Help

**For questions:**
- 💬 **GitHub Discussions** for general questions
- 🐛 **GitHub Issues** for bug reports
- 📧 **Email maintainers** for private matters

**For learning:**
- 📚 Read the documentation
- 🔍 Browse existing examples
- 🧪 Start with simple test cases

## 🎯 Development Priorities

### Current Focus Areas

**High Priority:**
1. **Validation and benchmarking**
2. **Documentation improvements**
3. **Platform compatibility**
4. **Performance optimization**

**Medium Priority:**
1. **New potential models**
2. **Additional examples**
3. **Build system improvements**

**Future Goals:**
1. **Parallelization** (MPI/OpenMP)
2. **Modern Fortran features**
3. **Python interface**

## 📚 Resources

### Physics Background
- Ichimura, Austern, Vincent formalism
- DWBA theory and applications
- Nuclear reaction theory

### Technical Resources
- [Fortran Standards](https://fortran-lang.org/)
- [Git Workflow](https://git-scm.com/book)
- [Scientific Computing Best Practices](https://software-carpentry.org/)

## 🙏 Recognition

Contributors will be:
- 📝 **Listed in CONTRIBUTORS.md**
- 🏷️ **Tagged in release notes**
- 💫 **Acknowledged in publications** (major contributions)

---

Thank you for contributing to SMOOTHIE! Your efforts help advance nuclear physics research worldwide. 🌟