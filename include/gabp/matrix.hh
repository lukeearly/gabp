#ifndef __MATRIX_HH__
#define __MATRIX_HH__

#include <memory>
#include <iostream>
#include <functional>
#include <algorithm>

/**
 * @brief The %gmat namespace includes the linear algebra backend for GaBP.
 */
namespace gmat {
    /**
    * @brief %matrix class for linear algebra behind inference algorithms.
    * @tparam T Type of elements.
    * @tparam m Number of rows.
    * @tparam n Number of columns.
    */
    template <typename T, size_t m, size_t n>
    class matrix {
    public:
        /**
        * @brief Gets the value of the element at coordinate i,j.
        * @param i Row coordinate.
        * @param j Column coordinate.
        * @return Value at i,j.
        *
        * @pre i < m
        * @pre j < n
        */
        virtual T get(size_t i, size_t j) const;

        /**
        * @brief Sets the value of the element at coordinate i,j.
        * @param i Row coordinate.
        * @param j Column coordinate.
        * @param value Value to be set.
        * @return New value.
        *
        * @pre i < m
        * @pre j < n
        */
        virtual T set(size_t i, size_t j, T value);

        /**
        * @brief Left %matrix multiplication.
        * @tparam o The width of the right %matrix.
        * @param right An @a n*o %matrix.
        * @return The @a m*o product of this @a m*n %matrix and the right @a n*o %matrix supplied.
        */
        template <size_t o>
        std::shared_ptr<matrix<T, m, o>> operator*(matrix<T, n, o>& right)
        {
            auto ret = std::make_shared<matrix<T, m, o>>();
            matmul(this, right, *ret);
            return ret;
        }

        /**
        * @brief Entryise %matrix summation.
        * @param right Another @a m*n %matrix.
        * @return The @a m*n sum of this @a m*n %matrix and the right @a m*n %matrix supplied.
        */
        std::shared_ptr<matrix<T, m, n>> operator+(matrix<T, m, n>& right)
        {
            auto ret = std::make_shared<matrix<T, m, n>>();
            matadd(this, right, *ret);
            return ret;
        }

        /**
         * @brief Compares this %matrix with another.
         * @param right Other %matrix.
         * @return true if each entry of the two matrices are equal.
         *         false if the two matrices differ.
         *
         * @warning This function checks strict equality, and is not recommended for
         *          floating point matrices. Use comppred with a thresholding predicate
         *          instead.
         */
        bool operator==(matrix<T, m, n>& right)
        {
            for (int i = 0; i < m; ++i) {
                for (int j = 0; j < n; ++j) {
                    if (this->get(i, j) != right.get(i, j)) {
                        return false;
                    }
                }
            }
            return true;
        }

        /**
         * @brief Compares this %matrix with another by a predicate.
         * @param right Other %matrix.
         * @param pred Predicate to compare two values of type T.
         * @return true if each entry of the two matrices pass the predicate.
         *         false if any entry of the two matrices fails the predicate.
         */
        bool cmppred(matrix<T, m, n>& right, std::function<bool(T, T)> pred)
        {
            for (int i = 0; i < m; ++i) {
                for (int j = 0; j < n; ++j) {
                    if (!pred(this->get(i, j), right.get(i, j))) {
                        return false;
                    }
                }
            }
            return true;
        }

        friend std::ostream& operator<<(std::ostream&out, const matrix<T, m, n>& mat)
        {
            out << "[ ";
            for (size_t j = 0; j < n; ++j) {
                out << "\t" << mat.get(0, j);
            }
            for (size_t i = 1; i < m; ++i) {
                out << "\n  ";
                for (size_t j = 0; j < n; ++j) {
                    out << "\t" << mat.get(i, j);
                }
            }
            out << "\t ]\n";
            return out;
        }

    private:
        /**
         * @brief Calculates the address of an element.
         * @param i 
         * @param j 
         * @return Pointer to element.
         */
        // virtual T* address(size_t i, size_t j);
    };

    template <typename T, size_t m, size_t n, size_t M, size_t N>
    class submatrix;

    template <typename T, size_t m, size_t n>
    class basematrix : public matrix<T, m, n> {
    public:
        /**
        * @brief Creates a %basematrix object.
        * @warning Not necessarily zero-valued.
        *
        * This constructor does not clear or set the array.
        */
        basematrix() { };

        /**
        * @brief Creates a %basematrix object with copies of an exemplar element.
        * @param ex Exemplar element.
        *
        * This constructor fills the %basematrix with @a m*n copies of ex.
        */
        basematrix(T ex)
        {
            std::fill_n((T*) m_elements, m * n, ex);
        };
        
        /**
        * @brief Creates a %basematrix object.
        * @param ptr Raw pointer to array of @a m*n elements in memory.
        * 
        * This constructor copies the values from ptr into the %basematrix.
        */
        basematrix(T* ptr)
        {
            std::copy_n(ptr, m * n, (T *) this->m_elements);
        };

        /**
        * @brief %basematrix copy constructor.
        * @param other Existing %basematrix of identical element type and dimensions.
        */
        basematrix(const basematrix<T, m ,n>& other)
        {
            std::copy_n(other.m_elements, m * n, this->m_elements);
        };

        /**
        * @brief %basematrix constructor from a submatrix.
        * @param other Existing %submatrix of identical element type and dimensions.
        * 
        * @todo Optimize this by copying rows or half-rows (need to 
        *       split if the submatrix wraps)
        */
        template<size_t M, size_t N>
        basematrix(const submatrix<T, m , n, M, N>& other)
        {
            //T* address = other->address(0, 0);
            for (int i = 0; i < m; ++i) {
                for (int j = 0; j < n; ++j) {
                    this->m_elements[i][j] = other.get(i, j);
                }
            }
        };

        /**
        * @brief Gets the value of the element at coordinate i,j.
        * @param i Row coordinate.
        * @param j Column coordinate.
        * @return Value at i,j.
        *
        * @pre i < m
        * @pre j < n
        */
        T get(size_t i, size_t j) const override
        {
            return m_elements[i][j];
        }

        /**
        * @brief Sets the value of the element at coordinate i,j.
        * @param i Row coordinate.
        * @param j Column coordinate.
        * @param value Value to be set.
        * @return New value.
        *
        * @pre i < m
        * @pre j < n
        */
        T set(size_t i, size_t j, T value) override
        {
            return m_elements[i][j] = value;
        }

        template<size_t sm, size_t sn>
        submatrix<T, sm, sn, m, n> submatrix(size_t i, size_t j);

    private:
        // T* address(size_t i, size_t j) override
        // {
        //     return &m_elements[i][j];
        // }

        /**
        * @brief 2-dimensional @a m*n array containing all elements of the %basematrix.
        */
        T m_elements[m][n];
    };

    /**
     * @brief Class that shadows any matrix object and represents a rectangular selection of it (wrapping at boundaries).
     * @tparam T Type of elements.
     * @tparam m Number of rows in %submatrix.
     * @tparam n Number of columns in %submatrix.
     * @tparam M Number of rows in original %matrix.
     * @tparam N Number of columns in original %matrix.
     */
    template <typename T, size_t m, size_t n, size_t M, size_t N>
    class submatrix : public matrix<T, m, n> {
    public:
        /**
         * @brief Creates a %submatrix object that directly mirrors a matrix.
         * @param mat %shared_ptr to the %matrix to be shadowed.
         */
        submatrix(std::shared_ptr<matrix<T, M, N>> mat) : m_parent(mat), m_i(0), m_j(0) { }

        /**
         * @brief Creates a %submatrix object that directly mirrors a matrix.
         * @param mat %shared_ptr to the %matrix to be shadowed.
         * @param i Vertical offset from top of %matrix.
         * @param j Horizontal offset from left of %matrix.
         */
        submatrix(std::shared_ptr<matrix<T, M, N>> mat, size_t i, size_t j) : m_parent(mat), m_i(i), m_j(j) { }

        T get(size_t i, size_t j) const override
        {
            return m_parent->get((i + m_i) % M, (j + m_j) % N);
        };

        T set(size_t i, size_t j, T value) override
        {
            return m_parent->set((i + m_i) % M, (j + m_j) % N, value);
        };

    private:
        // T* address(size_t i, size_t j) override
        // {
        //     return m_parent->address((i + m_i) % M, (j + m_j) % N);
        // }

        std::shared_ptr<matrix<T, M, N>> m_parent;
        size_t m_i, m_j;
    };

    /**
    * @brief Calculates the determinant of the %matrix.
    * @tparam T Type of elements.
    * @tparam n Number of rows and number of columns.
    * @param mat Shared pointer to %matrix to calculate determinant of.
    * @return Determinant of type T.
    */
    template <typename T, size_t n>
    T det(std::shared_ptr<matrix<T, n, n>> mat)
    {
        T acc = 0;
        bool plus = true;
        for (int j = 0; j < n; ++j) {
            std::shared_ptr<matrix<T, n - 1, n - 1>> sub = std::make_shared<submatrix<T, n - 1, n - 1, n, n>>(mat, 1, j + 1);
            auto sd = det(sub);
            if (plus) {
                acc += mat->get(0, j) * sd;
            } else {
                acc -= mat->get(0, j) * sd;
            }
            plus = !plus;
        }
        return acc;
    }

    template <typename T>
    T det(std::shared_ptr<matrix<T, 1, 1>>& mat)
    {
        return mat->get(0, 0);
    }

    /**
    * @brief Calculates the %inverse of a square %matrix and writes it into dest.
    * @tparam T Type of elements.
    * @tparam n Number of rows and number of columns.
    * @param src Reference to %matrix to invert.
    * @param dest Reference to %matrix to write results.
    * @return true if src is singular (ie non-invertible).
    *         false if src is non-singular (ie invertible).
    * 
    * @invariant src is unchanged.
    *
    * This function writes the inverse of src into dest, unless src is singular,
    * in which case it return true and leaves dest unchanged.
    */
    template <typename T, size_t n>
    bool inverse(matrix<T, n, n> &src, matrix<T, n, n> &dest)
    {
        return false;
    }

    /**
     * @brief Calculates the product of two matrices and writes it into dest.
     * @tparam T Type of elements.
     * @tparam m,n,o Dimensions of the three matrices involved.
     *               An @a m*n %matrix times an @a n*o %matrix is an @a m*o matrix.
     * @param left Reference to left %matrix to multiply.
     * @param right Reference to right %matrix to multiply.
     * @param dest Reference to %matrix to write results
     *
     * @invariant left, right are unchanged.
     */
    template <typename T, size_t m, size_t n, size_t o>
    void matmul(matrix<T, m, n>& left, matrix<T, n, o>& right, matrix<T, m, o>& dest)
    {
        for (int i = 0; i < m; ++i) {
            for (int j = 0; j < o; ++j) {
                T a = 0;
                for (int k = 0; k < n; ++k) {
                    a += left.get(i, k) * right.get(k, j);
                }
                dest.set(i, j, a);
            }
        }
    }

    /**
     * @brief Calculates the entrywise sum of two matrices and writes it to dest.
     * @tparam T Type of elements.
     * @tparam m,n Dimensions of the three matrices involved.
     * @param left Reference to left %matrix to add.
     * @param right Reference to right %matrix to add.
     * @param dest Reference to %matrix to write results
     *
     * @invariant left, right are unchanged.
     */
    template <typename T, size_t m, size_t n>
    void matadd(matrix<T, m, n>& left, matrix<T, m, n>& right, matrix<T, m, n>& dest)
    {
        for (int i = 0; i < m; ++i) {
            for (int j = 0; j < n; ++j) {
                dest.set(i, j, left.get(i, j) + right.get(i, j));
            }
        }
    }
}

#endif // __MATRIX_HH__
