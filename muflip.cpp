#include <GLFW/glfw3.h>
#include <algorithm>
#include <cmath>
#include <cstring>

#define FORCC   for (int j=0; j<N; ++j) for (int i=0; i<N; ++i)
#define FORCCIN for (int j=1; j<N-1; ++j) for (int i=1; i<N-1; ++i)
#define FORCORNERS(i,j) for (int jj=(j); jj<(j)+2; ++jj) for (int ii=(i); ii<(i)+2; ++ii)
#define IDXX(i,j) ((i)+(j)*Np1)
#define IDXY(i,j) (Y+(i)+(j)*N)

enum Cell { EMPTY, FLUID, SOLID };

const int N = 32;       // the size of the world is N*N
const int Np1 = N+1;
const int Y = N*(N+1);  // skip in MAC grids to the Y faces
double u[2*N*(N+1)];    // grid velocities
double ux[2*N*(N+1)];   // old grid velocities
double m[2*N*(N+1)];    // masses on the grid
double p[N*N];          // pressure cells
double divu[N*N];       // divergence of u, and the RHS of the Poisson eqn
Cell flags[N*N];        // solid/fluid/empty flags

const int PN = (N-2)/4*(N-4)*4;             // number of particles
double px[PN], py[PN], pvx[PN], pvy[PN];    // particle positions and velocities

double dt = 0.01;

void init() { // Dam Break
    int idx = 0;
    memset(pvx, 0, sizeof(pvx));
    memset(pvy, 0, sizeof(pvy));
    FORCCIN
        if (idx < PN && i-1 < (N-2)/4)
            FORCORNERS(0,0) px[idx] = i+0.25+ii*0.5, py[idx++] = j+0.25+jj*0.5; // add 4 particles to each initial cell
}

void particles2grid() {
    memset(u, 0, sizeof(u));
    memset(m, 0, sizeof(m));
    FORCC flags[i+j*N] = (i==0 || j==0 || i==N-1 || j==N-1) ? SOLID : EMPTY;
    for (int k=0; k<PN; ++k) {
        int i(px[k]), j(py[k]), fi(px[k]-0.5), fj(py[k]-0.5);
        flags[i+j*N] = FLUID;
        FORCORNERS(i,fj) {
            u[IDXX(ii,jj)] += pvx[k]*(1-fabs(ii-px[k]))*(1-fabs(jj+0.5-py[k]));
            m[IDXX(ii,jj)] += (1-fabs(ii-px[k]))*(1-fabs(jj+0.5-py[k]));
        }
        FORCORNERS(fi,j) {
            u[IDXY(ii,jj)] += pvy[k]*(1-fabs(ii+0.5-px[k]))*(1-fabs(jj-py[k]));
            m[IDXY(ii,jj)] += (1-fabs(ii+0.5-px[k]))*(1-fabs(jj-py[k]));
        }
    }
    for (int k=0; k<2*N*Np1; ++k)
        if (m[k]>1e-8) { m[k] = 1/m[k]; u[k] *= m[k]; }
    memcpy(ux, u, sizeof(u));
}

void pressureSolve() {
    FORCCIN { divu[i+j*N] = -(u[i+1+j*Np1]-u[i+j*Np1] + u[Y+i+(j+1)*N]-u[Y+i+j*N]); }
    for (int k=1; k<N-1; ++k) {
        divu[1+k*N] -= u[1+k*Np1];
        divu[N-2+k*N] += u[N-1+k*Np1];
        divu[k+N] -= u[Y+k+N];
        divu[k+(N-2)*N] += u[Y+k+(N-1)*N];
    }
    memset(p, 0, sizeof(p));
    for (int k=0; k<10*N*N; ++k)    // Gauss-Seidel solve for pressure
        FORCCIN {
            if (flags[i+j*N] != 1) continue;
            double l=i>1?-1:0, r=i<N-2?-1:0, b=j>1?-1:0, t=j<N-2?-1:0, c=-(l+r+b+t);
            p[i+j*N] = 1/c*(divu[i+j*N]-l*p[i-1+j*N]-r*p[i+1+j*N]-b*p[i+j*N-N]-t*p[i+j*N+N]);
        }
    for (int j=1; j<N-1; ++j)
        for (int i=1; i<N; ++i)
            u[i+j*Np1] -= p[i+j*N]-p[i-1+j*N];
    for (int j=1; j<N; ++j)
        for (int i=1; i<N-1; ++i)
            u[Y+i+j*N] -= p[i+j*N]-p[i+j*N-N];
    for (int k=1; k<N-1; ++k)
        u[1+k*Np1] = u[N-1+k*Np1] = u[Y+k+N] = u[Y+k+(N-1)*N] = 0;
}

void updateAndAdvect(double flip) {
    for (int k=0; k<PN; ++k) {
        int i(px[k]), j(py[k]), fi(px[k]-0.5), fj(py[k]-0.5);
        double vx(0), vy(0), vxOld(0), vyOld(0);
        FORCORNERS(i,fj) {
            vx += u[IDXX(ii,jj)]*(1-fabs(ii-px[k]))*(1-fabs(jj+0.5-py[k]));
            vxOld += ux[IDXX(ii,jj)]*(1-fabs(ii-px[k]))*(1-fabs(jj+0.5-py[k]));
        }
        FORCORNERS(fi,j) {
            vy += u[IDXY(ii,jj)]*(1-fabs(ii+0.5-px[k]))*(1-fabs(jj-py[k]));
            vyOld += ux[IDXY(ii,jj)]*(1-fabs(ii+0.5-px[k]))*(1-fabs(jj-py[k]));
        }
        pvx[k] = (1-flip)*vx + flip*(pvx[k]+vx-vxOld);
        pvy[k] = (1-flip)*vy + flip*(pvy[k]+vy-vyOld);
        px[k] = std::min(std::max(px[k]+vx*dt,1.001),N-1.001);
        py[k] = std::min(std::max(py[k]+vy*dt,1.001),N-1.001);
    }
}

int main() {
    if (!glfwInit()) exit(EXIT_FAILURE);
    GLFWwindow *window = glfwCreateWindow(800, 800, "muFLIP", NULL, NULL);
    glfwMakeContextCurrent(window);
    glPointSize(3);
    init();
    while (!glfwWindowShouldClose(window)) {
        particles2grid();
        for (int k=0; k<N*Np1; ++k)
            if (m[Y+k]>1e-8) u[Y+k] += -9.81*dt*(N-2);
        pressureSolve();
        updateAndAdvect(0.96);

        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
        glBegin(GL_POINTS);
        for (int k=0; k<PN; ++k)
            glVertex2d(1.9*(px[k]/N-0.5),1.9*(py[k]/N-0.5));
        glEnd();
        FORCC {
            if (flags[i+j*N] != SOLID) continue;
            glBegin(GL_LINE_LOOP);
            glVertex2d(1.9*(i/(double)N-0.5),1.9*(j/(double)N-0.5));
            glVertex2d(1.9*(i/(double)N-0.5),1.9*((j+1)/(double)N-0.5));
            glVertex2d(1.9*((i+1)/(double)N-0.5),1.9*((j+1)/(double)N-0.5));
            glVertex2d(1.9*((i+1)/(double)N-0.5),1.9*(j/(double)N-0.5));
            glEnd();
        }
        glfwSwapBuffers(window);
        glfwPollEvents();
    }
    return 0;
}
