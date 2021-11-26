#include "catch.hh"

#include <iostream>
#include <memory>
#include "gabp/matrix.hh"

TEST_CASE( "matrix instantiation" , " [matrix] ") {
    float m[3][3] = {
        {1.0, 0.0, 4.0},
        {0.0, 2.0, 0.0},
        {0.0, 0.0, 3.0}
    };
    auto mat = std::make_shared<gmat::basematrix<float, 3, 3>>((float*) m);
    REQUIRE( mat->get(0, 0) == 1.0 );
    REQUIRE( mat->get(0, 2) == 4.0 );
    REQUIRE( mat->get(2, 2) == 3.0 );
}

TEST_CASE( "matrix summation", " [matrix] ") {
    int a[3][3] = {
        {40,  2,   98},
        {36,  15,  52},
        {52,  34,  77}
    };
    int b[3][3] = {
        {37,  97,  77},
        {29,  3,   75},
        {92,  6,   14}
    };
    int c[3][3] = {
        {77,  99,  175},
        {65,  18,  127},
        {144, 40,  91}
    };
    gmat::basematrix<int, 3, 3> ma((int*) a);
    gmat::basematrix<int, 3, 3> mb((int*) b);
    gmat::basematrix<int, 3, 3> mc((int*) c);
    gmat::basematrix<int, 3, 3> msum;
    gmat::matadd(ma, mb, msum);
    auto good = ( msum == mc );
    REQUIRE( good );
}

TEST_CASE( "matrix multiplication", " [matrix] ") {
    int a[4][3] = {
        {1, 2, 3},
        {4, 5, 6},
        {7, 8 ,9},
        {10, 11, 12}
    };
    int b[3][2] = {
        {13, 14},
        {15, 16},
        {17, 18}
    };
    int c[4][2] = {
        {94, 100},
        {229, 244},
        {364, 388},
        {499, 532}
    };
    gmat::basematrix<int, 4, 3> ma((int*) a);
    gmat::basematrix<int, 3, 2> mb((int*) b);
    gmat::basematrix<int, 4, 2> mc((int*) c);
    gmat::basematrix<int, 4, 2> mprod;
    gmat::matmul(ma, mb, mprod);
    auto good = ( mprod == mc );
    REQUIRE( good );
}

TEST_CASE( "matrix determinant", "[matrix]" ) {
    SECTION( "1x1 matrix" ) {
        std::shared_ptr<gmat::matrix<int, 1, 1>> ma = std::make_shared<gmat::basematrix<int, 1, 1>>(2);
        REQUIRE( gmat::det(ma) == 2 );
    }

    SECTION( "2x2 matrix" ) {
        int b[2][2] = {
            {7, 13},
            {18, 6}
        };

        std::shared_ptr<gmat::matrix<int, 2, 2>> mb = std::make_shared<gmat::basematrix<int, 2, 2>>((int *)b);
        REQUIRE( gmat::det(mb) == 7 * 6 - 13 * 18 );
    }
}