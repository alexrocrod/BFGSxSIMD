#pragma once

#include <iterator>
#include <iostream>
#include <random>

#include <aadc/aadc.h>

extern "C" AADC_API void CAAD_MatrixVectorProduct(
    const double* mat, const double* vec, double* res,
    uint64_t num_r, uint64_t num_c
);

namespace aadc {

template<class mdouble>
class Matrix {
public:
    Matrix(uint64_t nr, uint64_t nc)
        : m_nr(nr), m_nc(nc)
        , m_data(m_nr * m_nc)
    {}
    Matrix(uint64_t nr, uint64_t nc, mdouble val)
        : m_nr(nr), m_nc(nc)
        , m_data(m_nr * m_nc, val)
    {}
    Matrix() : Matrix(0,0) {}
public:
    Matrix<mdouble>& operator = (const Matrix& other) {
        m_nr = other.m_nr;
        m_nc = other.m_nc;
        m_data = other.m_data;
        return *this;
    }

    mdouble& operator ()(const uint64_t r, const uint64_t c) {
        return m_data[r * m_nc + c];
    }
    const mdouble& operator ()(const uint64_t r, const uint64_t c) const {
        return m_data[r * m_nc + c];
    }

    mdouble* operator [](const uint64_t r) {
        return &m_data[r * m_nc];
    }
    const mdouble* operator [](const uint64_t r) const {
        return &m_data[r * m_nc];
    }

    void resize(uint64_t nr, uint64_t nc) {
        m_nr = nr; m_nc = nc;
        m_data.resize(nr*nc);
    }

    void resize(uint64_t nr, uint64_t nc, const mdouble& init_val) {
        m_nr = nr; m_nc = nc;
        m_data.resize(nr*nc, init_val);
    }
    
    uint64_t rows() const { return m_nr; }
    uint64_t cols() const { return m_nc; }

private:
    uint64_t m_nr, m_nc;
    std::vector<mdouble> m_data;
};

typedef Matrix<idouble> iMatrix;
typedef Matrix<double> ScalarMatrix;
template <class mmType>
using AVXMatrix = std::vector<mmVector<mmType> >;
typedef Matrix<aadc::AADCArgument> MatrixArg;
typedef Matrix<aadc::AADCScalarArgument> ScalarMatrixArg;
typedef Matrix<aadc::AADCResult> MatrixResult;

typedef std::vector<idouble> iVector;
typedef std::vector<double> ScalarVector;
template <class mmType>
using AVXVector = mmVector<mmType>;
typedef std::vector<aadc::AADCArgument> VectorArg;
typedef std::vector<aadc::AADCResult> VectorRes;
typedef std::vector<aadc::AADCScalarArgument> ScalarVectorArg;

template<class Matrix>
void initNNMatrix( Matrix& m, int size_row, int size_col, int seed = 0, double stddev = 1.) {
    std::mt19937_64 gen(seed);
    std::normal_distribution<> normal_dist(0, stddev);
    m.resize(size_row, size_col);
    for (uint64_t ri = 0; ri < size_row; ++ri) {
        for (uint64_t ci = 0; ci < size_col; ++ci) {
            m[ri][ci] = normal_dist(gen);
        }
    }
}

template<class Vector>
void initNNVector( Vector& m, uint64_t size) {
    m.resize(size);
    for (uint64_t ri = 0; ri < size; ++ri) m[ri] = 0.;
}

template<class mmType>
void initAVXMatrix( AVXMatrix<mmType>& m, uint64_t size_row, uint64_t size_col) {
    m.resize(size_row, size_col);
    for (uint64_t ri = 0; ri < size_row; ++ri) {
        for (uint64_t ci = 0;ci < size_col; ++ci) {
            m[ri][ci] = mmSetConst<mmType>(0.1);
        }
    }
}

inline iMatrix transpose(const iMatrix& m) {
    iMatrix matrix; 
    matrix.resize(m.rows(), m.cols());
    for (uint64_t ri = 0; ri < m.rows(); ++ri) {
        for (uint64_t ci = 0;ci < m.cols(); ++ci) {
            matrix[ri][ci] = m[ci][ri];
        }
    }
    return matrix;
}


template<class mmType>
void compareMatrixResults(const AVXMatrix<mmType>& mv, const ScalarMatrix& m) {
    if (mv.rows() != m.rows()) { std::cout << "m.rows,mv.rows = " << mv.rows() << " , " << m.rows() << "\n";}
    if (mv.cols() != m.cols()) { std::cout << "m.cols,mv.cols = " << mv.cols() << " , " << m.cols() << "\n";}
    for (uint64_t ri = 0; ri < mv.rows(); ++ri) {
        for (uint64_t ci = 0; ci < mv.cols(); ++ci) {
            if (abs(mmSum(mv[ri][ci])-m[ri][ci]) > 1e-10) {
                std::cout 
                    << "ri, ci, sum(mv), m = " << ri << " , " << ci << 
                    " , " << mmSum(mv[ri][ci]) << " , " << m[ri][ci] << "\n"
                ;
            }  
        }
    }
}

template<class mmType>
void compareVectorResults(const AVXVector<mmType>& vv, const ScalarVector& v) {
    if (vv.size() != v.size()) {std::cout << "m.size,mv.size = " << vv.size() << " , " << v.size() << "\n";}
    for (uint64_t ri = 0; ri < vv.size(); ++ri) {
        if (abs(mmSum(vv[ri])-v[ri]) >  1e-10) {
            std::cout << "ri,  sum(vv[ri]), v[ri] = " << ri << " , "  << mmSum(vv[ri]) << " , " << v[ri] << "\n";
        }
    }
}

template<class mmType, typename scalarType>  
inline void restructurizeData (
    const std::vector<std::vector<scalarType>>& scalar_data, 
    std::vector<AVXVector<mmType>>& mm_data
) {
    int AVX_size = aadc::mmSize<mmType>();
    mm_data.resize(scalar_data.size()/AVX_size);
    for (uint64_t i = 0; i < mm_data.size(); ++i) {
        mm_data[i].resize(scalar_data[0].size());
        for (uint64_t j=0; j<scalar_data[0].size(); j++) {
            for (uint64_t ci = 0; ci < AVX_size; ++ci) {
                toDblPtr(mm_data[i][j])[ci]=scalar_data[i*AVX_size+ci][j];
            }
        }
    }
}

inline void markMatrixAsInput(MatrixArg& mi, iMatrix& m, const bool diff) {
    uint64_t addr_start(CAAD_iReserveContinuousAddressSpace(false, m.rows()*m.cols()));
    mi.resize(m.rows(), m.cols());
    for (uint64_t ri = 0; ri < m.rows(); ++ri) {
        for (uint64_t ci = 0;ci < m.cols(); ++ci) {
            m[ri][ci].initVarAtAddress(addr_start);++addr_start;
            if (diff) {
                mi[ri][ci] = m[ri][ci].markAsInput();
            } else mi[ri][ci] = m[ri][ci].markAsInputNoDiff();
        }
    }
}


inline void markMatrixAsOutput(MatrixResult& mi, iMatrix& m) {
    mi.resize(m.rows(), m.cols());
    for (uint64_t ri = 0; ri < m.rows(); ++ri) {
        for (uint64_t ci = 0;ci < m.cols(); ++ci) {
            mi[ri][ci] = m[ri][ci].markAsOutput();
        }
    }
}


inline void markScalarMatrixAsInput(ScalarMatrixArg& mi, iMatrix& m) {
    uint64_t avx_size = CAAD_GetAVXSize();
    uint64_t rem(m.rows() % avx_size);
    uint64_t avx_r = (m.rows() / avx_size) + (rem ? 1 : 0);
    uint64_t nr = m.rows(), nc = m.cols();

    //uint64_t addr_start(CAAD_iReserveContinuousAddressSpace(true, m.size()*m[0].size()));
    uint64_t addr_start(CAAD_iReserveContinuousAddressSpace(true, avx_size*avx_r*nc));
    mi.resize(m.rows(), m.cols());
    for (uint64_t avx_ri = 0; avx_ri < avx_r; ++avx_ri) {
        for (uint64_t avx_i = 0; avx_i < avx_size; ++avx_i) {
            uint64_t ri = avx_ri * avx_size + avx_i;
            if (ri < nr) {
                for (uint64_t ci = 0;ci < nc; ++ci) {
//                    uint64_t addr(addr_start + avx_ri*avx_size*nc+avx_i+ci*avx_size);
                    m[ri][ci].initScalarVarAtAddress(addr_start + avx_ri*avx_size*nc+avx_i+ci*avx_size);
                    mi[ri][ci] = m[ri][ci].markAsScalarInput();
                }
            }
        }
    }
}

inline void markVectorAsInput(VectorArg& vi, const iVector& v, const bool Diff) {
    uint64_t addr_start(CAAD_iReserveContinuousAddressSpace(false, v.size()));
    vi.resize(v.size());
    for (uint64_t ri = 0; ri < v.size(); ++ri) {
            v[ri].initVarAtAddress(addr_start); ++addr_start;
            if (Diff) {
                vi[ri] = v[ri].markAsInput();
            } else {
                vi[ri] = v[ri].markAsInputNoDiff();
            }   
    }
}

inline void markVectorAsOutput(VectorRes& mi, const iVector& m) {
    mi.resize(m.size());
    for (uint64_t ri = 0; ri < m.size(); ++ri) {
        mi[ri] = m[ri].markAsOutput();
    }
}

inline void markVectorAsDiff(VectorArg& vi, const iVector& v) {
    uint64_t addr_start(CAAD_iReserveContinuousAddressSpace(false, v.size()));
    vi.resize(v.size());
    for (uint64_t ri = 0; ri < v.size(); ++ri) {
        v[ri].initVarAtAddress(addr_start); ++addr_start;
        vi[ri] = v[ri].markAsDiff();
    }
}

inline void markScalarVectorAsInput(ScalarVectorArg& vi, iVector& v) {
    uint64_t addr_start(CAAD_iReserveContinuousAddressSpace(true, v.size()));
    vi.resize(v.size());
    for (uint64_t ri = 0; ri < v.size(); ++ri) {
        v[ri].initScalarVarAtAddress(addr_start); ++addr_start;
        vi[ri] = v[ri].markAsScalarInput();
    }
}

template<class mmType>
void getAVXVectorGradient(aadc::AADCWorkSpace<mmType>& ws, const VectorArg& vi, AVXVector<mmType>& v) {
    v.resize(vi.size());
    for (uint64_t ri = 0; ri < v.size(); ++ri) {
        v[ri] = ws.diff(vi[ri]);
    }
}

template<class mmType>
void getScalarVectorGradient(aadc::AADCWorkSpace<mmType>& ws, const ScalarVectorArg& vi, ScalarVector& v) {
    v.resize(vi.size());
    for (uint64_t ri = 0; ri < v.size(); ++ri) {
        v[ri] = ws.diff(vi[ri]);
    }
}

template<class mmType>
void getAVXMatrixGradient(aadc::AADCWorkSpace<mmType>& ws, const MatrixArg& mi, AVXMatrix<mmType>& m) {
    m.resize(mi.rows(), mi.cols());
    for (uint64_t ri = 0; ri < m.size(); ++ri) {
        for (uint64_t ci = 0;ci < m[ri].size(); ++ci) {
            m[ri][ci] = ws.diff(mi[ri][ci]);
        }
    }
}

template<class mmType>
inline void getScalarMatrixGradient(aadc::AADCWorkSpace<mmType>& ws, const ScalarMatrixArg& mi, ScalarMatrix& m) {
    m.resize(mi.rows(), mi.cols());
    for (uint64_t ri = 0; ri < m.rows(); ++ri) {
        for (uint64_t ci = 0;ci < m.cols(); ++ci) {
            m[ri][ci] = ws.diff(mi[ri][ci]);
        }
    }
}

template<class mmType>
inline void getScalarMatrixGradient(aadc::AADCWorkSpace<mmType>& ws, const ScalarMatrixArg& mi, iMatrix& m) {
    m.resize(mi.rows(), mi.cols());
    for (int ri = 0; ri < m.rows(); ++ri) {
        for (int ci = 0;ci < m.cols(); ++ci) {
            m[ri][ci] = ws.diff(mi[ri][ci]);
        }
    }
}

template<class mmType>
inline void getScalarMatrix(aadc::AADCWorkSpace<mmType>& ws, const ScalarMatrixArg& mi, iMatrix& m) {
    for (int ri = 0; ri < m.rows(); ++ri) {
        for (int ci = 0; ci < m.cols(); ++ci) {
            m[ri][ci] = ws.val(mi[ri][ci]);
        }
    }
}

template<class mmType>
inline void getScalarVector(aadc::AADCWorkSpace<mmType>& ws, const ScalarVectorArg& mi, iVector& m) {
    for (int ri = 0; ri < m.size(); ++ri) m[ri] = ws.val(mi[ri]);
}


template<class mmType>
void setScalarMatrixInput(aadc::AADCWorkSpace<mmType>& ws, const ScalarMatrixArg& mi, const ScalarMatrix& m) {
    for (uint64_t ri = 0; ri < m.rows(); ++ri) {
        for (uint64_t ci = 0;ci < m.cols(); ++ci) {
            ws.setVal(mi[ri][ci], m[ri][ci]);
        }
    }
}

template<class mmType>
void setScalarMatrixDiffZero(aadc::AADCWorkSpace<mmType>& ws, const ScalarMatrixArg& mi) {
    for (uint64_t ri = 0; ri < mi.rows(); ++ri) {
        for (uint64_t ci = 0;ci < mi.cols(); ++ci) {
            ws.diff(mi[ri][ci]) = 0.;
        }
    }
}

template<class mmType>
void multScalarMatrixByCoeff(aadc::AADCWorkSpace<mmType>& ws, const ScalarMatrixArg& mi, const double coeff) {
    for (uint64_t ri = 0; ri < mi.rows(); ++ri) {
        for (uint64_t ci = 0;ci < mi.cols(); ++ci) {
            ws.diff(mi[ri][ci]) *= coeff;
        }
    }
}

template<class mmType>
void multScalarVectorByCoeff(aadc::AADCWorkSpace<mmType>& ws, const ScalarVectorArg& mi, const double coeff) {
    for (uint64_t ri = 0; ri < mi.size(); ++ri) ws.diff(mi[ri]) *= coeff;
}

template<class mmType>
inline void setScalarVectorInput(aadc::AADCWorkSpace<mmType>& ws, const ScalarVectorArg& vi, const ScalarVector& v) {
    for (uint64_t ri = 0; ri < v.size(); ++ri) {
        ws.setVal(vi[ri], v[ri]);
    }
}

template<class mmType>
void setScalarVectorDiffZero(aadc::AADCWorkSpace<mmType>& ws,  const ScalarVectorArg& vi) {
    for (uint64_t ri = 0; ri < vi.size(); ++ri) ws.diff(vi[ri]) = 0.;
}

template<class mmType>  
inline void setAVXVector(aadc::AADCWorkSpace<mmType>& ws, const VectorArg& vi, const AVXVector<mmType>& v) {
    for (uint64_t ri = 0; ri < v.size(); ++ri) {
        ws.setVal(vi[ri], v[ri]);
    }
}
template<class mmType>
inline void setScalarVectorToAVX(aadc::AADCWorkSpace<mmType>& ws, const VectorArg& vi, const ScalarVector& v) {
    for (uint64_t ri = 0; ri < v.size(); ++ri) {
        ws.setVal(vi[ri], mmSetConst<mmType>(v[ri]));
    }
}

template<class mmType>
inline void setScalarMatrixInput(aadc::AADCWorkSpace<mmType>& ws, const ScalarMatrixArg& mi, const iMatrix& m) {
    for (int ri = 0; ri < m.rows(); ++ri) {
        for (int ci = 0;ci < m.cols(); ++ci) {
            ws.setVal(mi[ri][ci], m[ri][ci].val);
        }
    }
}


template<class mmType>
inline void setAVXMatrixInput(aadc::AADCWorkSpace<mmType>& ws, const MatrixArg& mi, const AVXMatrix<mmType> m) {
    for (uint64_t ri = 0; ri < m.rows(); ++ri) {
        for (uint64_t ci = 0;ci < m.cols(); ++ci) {
            ws.setVal(mi[ri][ci], m[ri][ci]);
        }
    }
}

template<class mmType>
inline void setScalarMatrixInputToAVX(aadc::AADCWorkSpace<mmType>& ws, const MatrixArg& mi, const ScalarMatrix& m) {
    for (uint64_t ri = 0; ri < m.rows(); ++ri) {
        for (uint64_t ci = 0;ci < m.cols(); ++ci) {
            ws.setVal(mi[ri][ci], mmSetConst<mmType>(m[ri][ci]));
        }
    }
}

// Experimental version. Works and optimized for situation when xs is markedAsScalar
inline void matrixVectorMult( const iMatrix& xs, const iVector& x, iVector& y ) {
    bool saved = idouble::recording;
    if (idouble::recording) {
        CAAD_MatrixVectorProduct(&(xs[0][0].val), &(x[0].val), &(y[0].val), xs.rows(), xs.cols());
    }
    idouble::recording = false;
    for (int i = 0; i < xs.rows(); i++) {
        idouble s = 0.0;
        for (int k = 0; k < xs.cols(); k++) s += xs[i][k] * x[k];
        y[i] = s;
    }
    idouble::recording = saved;
    return;
}

inline void matrixVectorMult(const Matrix<double>& xs, const ScalarVector& x, ScalarVector& y) {
    for (int i = 0; i < xs.rows(); i++) {
        double s = 0.;
        for (int k = 0; k < xs.cols(); k++) s += xs[i][k] * x[k];
        y[i] = s;
    }
    return;
}

inline std::ostream& operator<<( std::ostream& out, const iVector& v) {
  std::copy( v.begin(), v.end(), std::ostream_iterator<idouble>( out, " " ) );
  return out;
}

inline std::ostream& operator<<( std::ostream& out, const ScalarVector& v) {
  std::copy( v.begin(), v.end(), std::ostream_iterator<double>( out, " " ) );
  return out;
}

} // namespace aadc
