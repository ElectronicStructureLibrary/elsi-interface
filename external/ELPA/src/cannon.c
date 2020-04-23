#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <mpi.h>

#ifdef Add_
#define numroc numroc_
#define dlacpy dlacpy_
#define pdtran pdtran_
#define dlacpy dlacpy_
#define pdlacpy pdlacpy_
#define slacpy slacpy_
#define pstran pstran_
#define pslacpy pslacpy_
#define zlacpy zlacpy_
#define pztranc pztranc_
#define pzlacpy pzlacpy_
#define clacpy clacpy_
#define pctranc pctranc_
#define pclacpy pclacpy_
#endif

int numroc_(int*, int*, int*, int*, int*);

void pdlacpy_(char*, int*, int*, double*, int*, int*, int*, double*, int*, int*, int*);
void pdtran_(int*, int*, double*, double*, int*, int*, int*, double*, double*, int*, int*, int*);

void pslacpy_(char*, int*, int*, float*, int*, int*, int*, float*, int*, int*, int*);
void pstran_(int*, int*, float*, float*, int*, int*, int*, float*, float*, int*, int*, int*);

void pzlacpy_(char*, int*, int*, double complex*, int*, int*, int*, double complex*, int*, int*, int*);
void pztranc_(int*, int*, double complex*, double complex*, int*, int*, int*, double complex*, double complex*, int*, int*, int*);

void pclacpy_(char*, int*, int*, float complex*, int*, int*, int*, float complex*, int*, int*, int*);
void pctranc_(int*, int*, float complex*, float complex*, int*, int*, int*, float complex*, float complex*, int*, int*, int*);

void dlacpy(char*, int*, int*, double*, int*, double*, int*);
void dgemm_(char*, char*, int*, int*, int*, double*, double*, int*, double*, int*, double*, double*, int*);

void slacpy(char*, int*, int*, float*, int*, float*, int*);
void sgemm_(char*, char*, int*, int*, int*, float*, float*, int*, float*, int*, float*, float*, int*);

void zlacpy(char*, int*, int*, double complex*, int*, double complex*, int*);
void zgemm_(char*, char*, int*, int*, int*, double complex*, double complex*, int*, double complex*, int*, double complex*, double complex*, int*);

void clacpy(char*, int*, int*, float complex*, int*, float complex*, int*);
void cgemm_(char*, char*, int*, int*, int*, float complex*, float complex*, int*, float complex*, int*, float complex*, float complex*, int*);
int numroc(int*, int*, int*, int*, int*);

void pdlacpy(char*, int*, int*, double*, int*, int*, int*, double*, int*, int*, int*);
void pdtran(int*, int*, double*, double*, int*, int*, int*, double*, double*, int*, int*, int*);

void pslacpy(char*, int*, int*, float*, int*, int*, int*, float*, int*, int*, int*);
void pstran(int*, int*, float*, float*, int*, int*, int*, float*, float*, int*, int*, int*);

void pzlacpy(char*, int*, int*, double complex*, int*, int*, int*, double complex*, int*, int*, int*);
void pztranc(int*, int*, double complex*, double complex*, int*, int*, int*, double complex*, double complex*, int*, int*, int*);

void pclacpy(char*, int*, int*, float complex*, int*, int*, int*, float complex*, int*, int*, int*);
void pctranc(int*, int*, float complex*, float complex*, int*, int*, int*, float complex*, float complex*, int*, int*, int*);

void cannons_reduction_d(double* A, double* U, int np_rows, int np_cols, int my_prow, int my_pcol,
                         int* a_desc, double *Res, int ToStore, MPI_Comm row_comm, MPI_Comm col_comm)
{
   int na, nblk, i, j, Size_send_A, Size_receive_A, Size_send_U, Size_receive_U, Buf_rows, Buf_cols, where_to_send_A, from_where_to_receive_A, where_to_send_U, from_where_to_receive_U, last_proc_row, last_proc_col, cols_in_buffer_A, rows_in_buffer_A, intNumber;
   double *Buf_to_send_A, *Buf_to_receive_A, *Buf_to_send_U, *Buf_to_receive_U, *data_ptr, *Buf_A, *Buf_pos, *U_local_start, *Res_ptr, *M, *M_T, *A_local_start, *U_local_start_curr, *U_stored, *CopyTo, *CopyFrom, *U_to_calc;
   int ratio, num_of_iters, cols_in_buffer, rows_in_block, rows_in_buffer, curr_col_loc, cols_in_block, curr_col_glob, curr_row_loc, Size_receive_A_now, Nb, owner, cols_in_buffer_A_now;
   int Size_receive_A_nowMPI, Size_receive_AMPI, Size_receive_UMPI;

   int row_of_origin_U, rows_in_block_U, num_of_blocks_in_U_buffer, k, startPos, cols_in_buffer_U, rows_in_buffer_U, col_of_origin_A, curr_row_loc_res, curr_row_loc_A, curr_col_glob_res;
   int curr_col_loc_res, curr_col_loc_buf, proc_row_curr, curr_col_loc_U, A_local_index, LDA_A, LDA_A_new, index_row_A_for_LDA, ii, rows_in_block_U_curr, width, row_origin_U, rows_in_block_A, cols_in_buffer_A_my_initial, rows_in_buffer_A_my_initial, proc_col_min;
   int *SizesU;
   int Size_U_skewed, Size_U_stored, Curr_pos_in_U_stored, rows_in_buffer_A_now;
   double done = 1.0;
   double dzero = 0.0;
   int one = 1;
   int zero = 0;
   int na_rows, na_cols;

   MPI_Status status;
   MPI_Request request_A_Recv;
   MPI_Request request_A_Send;
   MPI_Request request_U_Recv;
   MPI_Request request_U_Send;

   na = a_desc[2];
   nblk = a_desc[4];
   na_rows = numroc(&na, &nblk, &my_prow, &zero, &np_rows);
   na_cols = numroc(&na, &nblk, &my_pcol, &zero, &np_cols);

   if(ToStore > (np_rows -1))
      if((my_prow == 0)&&(my_pcol == 0))
         printf("Buffering level is larger than (np_rows-1) !!!\n");
   if((my_prow == 0)&&(my_pcol == 0))
         printf("Buffering level = %d\n", ToStore);

   if (np_cols%np_rows != 0)
   {
      if((my_prow == 0)&& (my_pcol ==0))
         printf("!!!!! np_cols must be a multiple of np_rows!!!!! I do nothing! \n");
      return;
   }
   if (np_cols < np_rows != 0)
   {
      if((my_prow == 0)&& (my_pcol ==0))
         printf("np_cols < np_rows \n");
      return;
   }

   ratio = np_cols/np_rows;
   last_proc_row = ((na-1)/nblk) % np_rows;
   last_proc_col = ((na-1)/nblk) % np_cols;

   if(na%nblk == 0)
      if(my_pcol <= last_proc_col)
         Buf_cols = na_cols;
      else
         Buf_cols = na_cols + nblk;
   else
      if(my_pcol < last_proc_col)
         Buf_cols = na_cols;
      else if(my_pcol > last_proc_col)
         Buf_cols = na_cols + nblk;
      else
         Buf_cols = na_cols + nblk - na_cols%nblk;

  if(na%nblk == 0)
      if(my_prow <= last_proc_row)
         Buf_rows = na_rows;
      else
         Buf_rows = na_rows + nblk;
   else
      if(my_prow < last_proc_row)
         Buf_rows = na_rows;
      else if(my_prow > last_proc_row)
         Buf_rows = na_rows + nblk;
      else
         Buf_rows = na_rows + nblk - na_rows%nblk;

   intNumber = ceil((double)na/(double)(np_cols*nblk));
   Size_U_stored = ratio*nblk*nblk*intNumber*(intNumber+1)/2 + 2;

   U_stored = malloc((Size_U_stored*(ToStore+1))*sizeof(double));
   SizesU = malloc(ToStore*sizeof(int));
   Buf_to_send_A = malloc(ratio*Buf_cols*Buf_rows*sizeof(double));
   Buf_to_receive_A = malloc(ratio*Buf_cols*Buf_rows*sizeof(double));
   Buf_to_send_U = malloc(Size_U_stored*sizeof(double));
   Buf_to_receive_U = malloc(Size_U_stored*sizeof(double));
   if(ratio != 1)
      Buf_A = malloc(Buf_cols*Buf_rows*sizeof(double));
   M = malloc(na_rows*na_cols*sizeof(double));
   M_T = malloc(na_rows*na_cols*sizeof(double));
   for(i = 0; i < na_rows*na_cols; i++)
      M[i] = 0;

   if(ratio != 1)
      dlacpy("A", &na_rows, &na_cols, A, &na_rows, Buf_to_send_A, &na_rows);
   Size_receive_A = 0;

   for(i = 0; i < ratio; i++)
   {
      where_to_send_A = (my_pcol - my_prow - i*np_rows + np_cols)%np_cols;
      from_where_to_receive_A = (my_pcol + my_prow + i*np_rows)%np_cols;

      if(ratio != 1)
      {
         if(where_to_send_A != my_pcol)
         {
           MPI_Sendrecv(Buf_to_send_A, (int) (na_cols*na_rows), MPI_DOUBLE,(int) where_to_send_A, (int) zero, Buf_A, (int) (na_rows*Buf_cols), MPI_DOUBLE, (int) from_where_to_receive_A, (int) zero, row_comm, &status);
           MPI_Get_count(&status, MPI_DOUBLE, &Size_receive_A_nowMPI);
           Size_receive_A_now = (int) Size_receive_A_nowMPI;
           Size_receive_A_now = Size_receive_A_now/na_rows;
         }
         else
            Size_receive_A_now = na_cols;
         Size_receive_A = Size_receive_A + Size_receive_A_now;

         intNumber = from_where_to_receive_A/np_rows;

         CopyTo = &Buf_to_receive_A[intNumber*na_rows*nblk];
         if(where_to_send_A != my_pcol)
            CopyFrom = Buf_A;
         else
            CopyFrom = A;

         intNumber = ceil((double)Size_receive_A_now/(double)nblk);
         for(j = 0; j < intNumber; j++)
         {
            width = nblk;
            if(nblk*(j+1) > Size_receive_A_now)
               width = Size_receive_A_now - nblk*j;
            dlacpy("A", &na_rows, &width, CopyFrom, &na_rows, CopyTo, &na_rows);
            CopyTo = CopyTo + na_rows*nblk*ratio;
            CopyFrom = CopyFrom + na_rows*nblk;
         }
      }
      else
         if(my_prow > 0)
         {
            dlacpy("A", &na_rows, &na_cols, A, &na_rows, Buf_to_send_A, &na_rows);
            MPI_Sendrecv(Buf_to_send_A, (int) (na_cols*na_rows), MPI_DOUBLE, (int) where_to_send_A, (int) zero, Buf_to_receive_A, (int) (na_rows*Buf_cols), MPI_DOUBLE, (int) from_where_to_receive_A, (int) zero, row_comm, &status);
            MPI_Get_count(&status, MPI_DOUBLE, &Size_receive_AMPI);
            Size_receive_A = (int) Size_receive_AMPI;
     Size_receive_A = Size_receive_A/na_rows;
         }
         else
         {
            dlacpy("A", &na_rows, &na_cols, A, &na_rows, Buf_to_receive_A, &na_rows);
            Size_receive_A = na_cols;
         }
   }

   num_of_iters = ceil((double)na_cols/(double)nblk);

   where_to_send_U = (my_prow - my_pcol + np_cols)%np_rows;
   from_where_to_receive_U = (my_pcol + my_prow)%np_rows;

   if(where_to_send_U == my_prow)
      Buf_pos = Buf_to_receive_U;
   else
      Buf_pos = Buf_to_send_U;

   if(my_pcol >= my_prow)
      curr_col_loc = 0;
   else
      curr_col_loc = 1;

   num_of_iters = num_of_iters - curr_col_loc;
   curr_col_loc = curr_col_loc*nblk;

   if(my_pcol >= my_prow )
      rows_in_block = ceil(((double)(my_pcol + 1) - (double)my_prow)/(double)np_rows)*nblk;
   else
      rows_in_block = ratio*nblk;

   Size_send_U = 0;
   for(i = 0; i < num_of_iters; i++)
   {
      if(rows_in_block > na_rows)
         rows_in_block = na_rows;

      if ((na_cols - curr_col_loc) < nblk)
         cols_in_block = na_cols - curr_col_loc;
      else
         cols_in_block = nblk;

      if((rows_in_block > 0)&&(cols_in_block > 0))
      {
         data_ptr = &U[curr_col_loc*na_rows];
         dlacpy("A", &rows_in_block, &cols_in_block, data_ptr, &na_rows, Buf_pos, &rows_in_block);
         Buf_pos = Buf_pos + rows_in_block*cols_in_block;
         Size_send_U = Size_send_U + rows_in_block*cols_in_block;
      }
      curr_col_loc = curr_col_loc + nblk;
      rows_in_block = rows_in_block + ratio*nblk;
   }
   rows_in_buffer = rows_in_block - ratio*nblk;
   *Buf_pos = (double)rows_in_buffer;
   Size_send_U = Size_send_U + 1;

   if(where_to_send_U != my_prow)
   {

      MPI_Sendrecv(Buf_to_send_U, (int) Size_send_U, MPI_DOUBLE, (int) where_to_send_U, (int) zero, Buf_to_receive_U, (int) (Buf_rows*na_cols), MPI_DOUBLE, (int) from_where_to_receive_U, (int) zero, col_comm, &status);
      MPI_Get_count(&status, MPI_DOUBLE, &Size_receive_UMPI);
      Size_receive_U = (int) Size_receive_UMPI;
   }
   else
      Size_receive_U = Size_send_U;

   for(i = 0; i < Size_receive_U; i++)
      U_stored[i] = Buf_to_receive_U[i];
   Size_U_skewed = Size_receive_U;
   Curr_pos_in_U_stored = Size_U_skewed;

   where_to_send_A = (my_pcol - 1 + np_cols)%np_cols;
   from_where_to_receive_A = (my_pcol + 1)%np_cols;
   where_to_send_U = (my_prow - 1 + np_rows)%np_rows;
   from_where_to_receive_U = (my_prow + 1)%np_rows;

   for(j = 1; j < np_rows; j++)
   {

      data_ptr = Buf_to_send_A;
      Buf_to_send_A = Buf_to_receive_A;
      Buf_to_receive_A = data_ptr;

      data_ptr = Buf_to_send_U;
      Buf_to_send_U = Buf_to_receive_U;
      Buf_to_receive_U = data_ptr;

      Size_send_A = Size_receive_A;
      MPI_Isend(Buf_to_send_A, (int) (Size_send_A*na_rows), MPI_DOUBLE, (int) where_to_send_A, (int) zero, row_comm, &request_A_Send);
      MPI_Irecv(Buf_to_receive_A, (int) (Buf_cols*na_rows*ratio), MPI_DOUBLE, (int) from_where_to_receive_A, (int) zero, row_comm, &request_A_Recv);

      Size_send_U = Size_receive_U;
      MPI_Isend(Buf_to_send_U, (int) Size_send_U, MPI_DOUBLE, (int) where_to_send_U, (int) zero, col_comm, &request_U_Send);
      MPI_Irecv(Buf_to_receive_U, (int) (Buf_rows*na_cols), MPI_DOUBLE, (int) from_where_to_receive_U, (int) zero, col_comm, &request_U_Recv);

      rows_in_buffer = (int)Buf_to_send_U[Size_receive_U-1];
      row_origin_U = (my_pcol + my_prow + np_cols + j - 1)%np_rows;

      if((my_pcol >= my_prow)&&(my_pcol >= row_origin_U))
      {
         cols_in_buffer = na_cols;
         curr_col_loc_res = 0;
         curr_col_loc_buf = 0;
      }
      if((my_pcol < my_prow)&&(my_pcol < row_origin_U))
      {
         cols_in_buffer = na_cols - nblk;
         curr_col_loc_res = nblk;
         curr_col_loc_buf = 0;
      }
      if((my_pcol >= my_prow)&&(my_pcol < row_origin_U))
      {
         cols_in_buffer = na_cols - nblk;
         curr_col_loc_res = nblk;
         curr_col_loc_buf = 0;
      }
      if((my_pcol < my_prow)&&(my_pcol >= row_origin_U))
      {
         cols_in_buffer = na_cols;
         curr_col_loc_res = nblk;
         curr_col_loc_buf = nblk;
      }

      num_of_blocks_in_U_buffer = ceil(((double)cols_in_buffer - (double)curr_col_loc_buf)/(double)nblk);

      startPos = (curr_col_loc_buf + nblk)*curr_col_loc_buf/2;
      U_local_start = &Buf_to_send_U[startPos];
      Res_ptr = &M[curr_col_loc_res*na_rows];

      for (i = 0; i < num_of_blocks_in_U_buffer; i++)
      {
         curr_col_glob = (curr_col_loc_res/nblk)*nblk*np_cols + my_pcol*nblk;
         proc_row_curr = (curr_col_glob/nblk)%np_rows;
         rows_in_block_A = (curr_col_glob/(nblk*np_rows))*nblk;
         if(my_prow <= proc_row_curr)
            rows_in_block_A = rows_in_block_A + nblk;

         if(rows_in_block_A > na_rows)
            rows_in_block_A = na_rows;

         if((curr_col_loc_buf + nblk) <= cols_in_buffer)
            cols_in_block = nblk;
         else
            cols_in_block = cols_in_buffer - curr_col_loc_buf;

         rows_in_block_U = (curr_col_glob/(nblk*np_rows))*nblk;
         if(proc_row_curr >= row_origin_U)
            rows_in_block_U = rows_in_block_U + nblk;

         if(rows_in_block_U > rows_in_buffer)
            rows_in_block_U = rows_in_buffer;

         if ((rows_in_block_A > 0)&&(cols_in_block > 0))
            if (j == 1) {
               dgemm_("N", "N", &rows_in_block_A, &cols_in_block, &rows_in_block_U, &done, Buf_to_send_A, &na_rows, U_local_start, &rows_in_block_U, &dzero, Res_ptr, &na_rows);
     }
            else {
               dgemm_("N", "N", &rows_in_block_A, &cols_in_block, &rows_in_block_U, &done, Buf_to_send_A, &na_rows, U_local_start, &rows_in_block_U, &done, Res_ptr, &na_rows);
     }

         U_local_start = U_local_start + rows_in_block_U*cols_in_block;
         curr_col_loc_res = curr_col_loc_res + nblk;
         Res_ptr = &M[curr_col_loc_res*na_rows];
         curr_col_loc_buf = curr_col_loc_buf + nblk;
      }

      MPI_Wait(&request_A_Send, &status);
      MPI_Wait(&request_A_Recv, &status);

      MPI_Get_count(&status, MPI_DOUBLE, &Size_receive_AMPI);
      Size_receive_A = (int) Size_receive_AMPI;
      Size_receive_A = Size_receive_A / na_rows;

      MPI_Wait(&request_U_Send, &status);
      MPI_Wait(&request_U_Recv, &status);
      MPI_Get_count(&status, MPI_DOUBLE, &Size_receive_UMPI);
      Size_receive_U = (int) Size_receive_UMPI;

      if(j <= ToStore)
      {
         for(k = 0; k < Size_receive_U; k++)
            U_stored[Curr_pos_in_U_stored + k] = Buf_to_receive_U[k];
         Curr_pos_in_U_stored = Curr_pos_in_U_stored + Size_receive_U;
         SizesU[j-1] = Size_receive_U;
      }
   }

   rows_in_buffer = (int)Buf_to_receive_U[Size_receive_U-1];
   row_origin_U = (my_pcol + my_prow + np_cols + np_rows -1)%np_rows;

   if((my_pcol >= my_prow)&&(my_pcol >= row_origin_U))
   {
      cols_in_buffer = na_cols;
      curr_col_loc_res = 0;
      curr_col_loc_buf = 0;
   }
   if((my_pcol < my_prow)&&(my_pcol < row_origin_U))
   {
      cols_in_buffer = na_cols - nblk;
      curr_col_loc_res = nblk;
      curr_col_loc_buf = 0;
   }
   if((my_pcol >= my_prow)&&(my_pcol < row_origin_U))
   {
      cols_in_buffer = na_cols - nblk;
      curr_col_loc_res = nblk;
      curr_col_loc_buf = 0;
   }
   if((my_pcol < my_prow)&&(my_pcol >= row_origin_U))
   {
      cols_in_buffer = na_cols;
      curr_col_loc_res = nblk;
      curr_col_loc_buf = nblk;
   }

   num_of_blocks_in_U_buffer = ceil(((double)cols_in_buffer - (double)curr_col_loc_buf)/(double)nblk);

   startPos = (curr_col_loc_buf + nblk)*curr_col_loc_buf/2;
   U_local_start = &Buf_to_receive_U[startPos];
   Res_ptr = &M[curr_col_loc_res*na_rows];

   for (i = 0; i < num_of_blocks_in_U_buffer; i++)
   {
      curr_col_glob = (curr_col_loc_res/nblk)*nblk*np_cols + my_pcol*nblk;
      proc_row_curr = (curr_col_glob/nblk)%np_rows;
      rows_in_block_A = (curr_col_glob/(nblk*np_rows))*nblk;
      if(my_prow <= proc_row_curr)
         rows_in_block_A = rows_in_block_A + nblk;

      if(rows_in_block_A > na_rows)
         rows_in_block_A = na_rows;

      if((curr_col_loc_buf + nblk) <= cols_in_buffer)
         cols_in_block = nblk;
      else
         cols_in_block = cols_in_buffer - curr_col_loc_buf;

      rows_in_block_U = (curr_col_glob/(nblk*np_rows))*nblk;
      if(proc_row_curr >= row_origin_U)
         rows_in_block_U = rows_in_block_U + nblk;

      if(rows_in_block_U > rows_in_buffer)
         rows_in_block_U = rows_in_buffer;

      if ((rows_in_block_A > 0)&&(cols_in_block > 0))
         if (j == 1) {
            dgemm_("N", "N", &rows_in_block_A, &cols_in_block, &rows_in_block_U, &done, Buf_to_receive_A, &na_rows, U_local_start, &rows_in_block_U, &dzero, Res_ptr, &na_rows);
  }
         else {
            dgemm_("N", "N", &rows_in_block_A, &cols_in_block, &rows_in_block_U, &done, Buf_to_receive_A, &na_rows, U_local_start, &rows_in_block_U, &done, Res_ptr, &na_rows);
         }

      U_local_start = U_local_start + rows_in_block_U*cols_in_block;
      curr_col_loc_res = curr_col_loc_res + nblk;
      Res_ptr = &M[curr_col_loc_res*na_rows];
      curr_col_loc_buf = curr_col_loc_buf + nblk;
   }

   pdtran(&na, &na, &done, M, &one, &one, a_desc, &dzero, M_T, &one, &one, a_desc);

   if((ratio != 1)||(my_prow != 0))
      Buf_pos = Buf_to_send_A;
   else
      Buf_pos = Buf_to_receive_A;

   num_of_iters = ceil((double)na_cols/(double)nblk);

   cols_in_buffer_A_my_initial = 0;
   Size_send_A = 0;

   if(my_pcol <= my_prow)
   {
      curr_row_loc = 0;
      rows_in_buffer_A_my_initial = na_rows;
   }
   else
   {
      curr_row_loc = ceil((double)(((double)my_pcol - (double)my_prow)/(double)np_rows))*nblk;
      rows_in_buffer_A_my_initial = na_rows - curr_row_loc;
   }

   for(i = 0; i < num_of_iters; i++)
   {
      curr_col_loc = i*nblk;
      rows_in_block = na_rows - curr_row_loc;

      if ((na_cols - curr_col_loc) < nblk)
         cols_in_block = na_cols - curr_col_loc;
      else
         cols_in_block = nblk;

      if((rows_in_block > 0)&&(cols_in_block > 0))
      {
         A_local_start = &M_T[curr_col_loc*na_rows + curr_row_loc];
         dlacpy("A", &rows_in_block, &cols_in_block, A_local_start, &na_rows, Buf_pos, &rows_in_block);
         Buf_pos = Buf_pos + rows_in_block*cols_in_block;
         Size_send_A = Size_send_A + rows_in_block*cols_in_block;
         cols_in_buffer_A_my_initial = cols_in_buffer_A_my_initial + cols_in_block;
      }
      curr_row_loc = curr_row_loc + ratio*nblk;
   }
   *Buf_pos = (double)cols_in_buffer_A_my_initial;
   Size_send_A = Size_send_A + 1;

   proc_col_min = np_cols;
   for(i = 0; i < ratio; i++)
   {
      from_where_to_receive_A = (my_pcol + my_prow + i*np_rows)%np_cols;
      if(from_where_to_receive_A < proc_col_min)
         proc_col_min = from_where_to_receive_A;
   }

   Size_receive_A = 0;
   cols_in_buffer_A = 0;
   rows_in_buffer_A = 0;
   for(i = 0; i < ratio; i++)
   {
      where_to_send_A = (my_pcol - my_prow - i*np_rows + np_cols)%np_cols;
      from_where_to_receive_A = (my_pcol + my_prow + i*np_rows)%np_cols;

      if(ratio != 1)
      {
         if(where_to_send_A != my_pcol)
         {
            MPI_Sendrecv(Buf_to_send_A, (int) Size_send_A, MPI_DOUBLE, (int) where_to_send_A, (int) zero, Buf_A, (int) Size_U_stored, MPI_DOUBLE, (int) from_where_to_receive_A, (int) zero, row_comm, &status);
            MPI_Get_count(&status, MPI_DOUBLE, &Size_receive_A_nowMPI);
            Size_receive_A_now = (int) Size_receive_A_nowMPI;

            Size_receive_A = Size_receive_A + Size_receive_A_now - 1;

            cols_in_buffer_A_now = Buf_A[Size_receive_A_now-1];
            cols_in_buffer_A = cols_in_buffer_A + cols_in_buffer_A_now;

            if(from_where_to_receive_A <= my_prow)
            {
               rows_in_buffer_A_now = na_rows;
            }
            else
            {
               rows_in_buffer_A_now = na_rows - ceil((double)(((double)from_where_to_receive_A - (double)my_prow)/(double)np_rows))*nblk;
            }
            if(rows_in_buffer_A < rows_in_buffer_A_now)
               rows_in_buffer_A = rows_in_buffer_A_now;

            intNumber = from_where_to_receive_A/np_rows;
            if(proc_col_min <= my_prow)
               CopyTo = &Buf_to_receive_A[nblk*(na_rows*intNumber - nblk*(intNumber-1)*intNumber/2)];
            else
               CopyTo = &Buf_to_receive_A[nblk*(na_rows*intNumber - nblk*intNumber*(intNumber+1)/2)];
            CopyFrom = Buf_A;
         }
         else
         {
            cols_in_buffer_A_now = cols_in_buffer_A_my_initial;
            cols_in_buffer_A = cols_in_buffer_A + cols_in_buffer_A_now;

            rows_in_buffer_A_now = rows_in_buffer_A_my_initial;
            if(rows_in_buffer_A < rows_in_buffer_A_now)
               rows_in_buffer_A = rows_in_buffer_A_now;

            intNumber = my_pcol/np_rows;
            if(proc_col_min <= my_prow)
               CopyTo = &Buf_to_receive_A[nblk*(na_rows*intNumber - nblk*(intNumber-1)*intNumber/2)];
            else
               CopyTo = &Buf_to_receive_A[nblk*(na_rows*intNumber - nblk*intNumber*(intNumber+1)/2)];
            CopyFrom = Buf_to_send_A;

            Size_receive_A = Size_receive_A + Size_send_A - 1;
         }

         intNumber = ceil((double)cols_in_buffer_A_now/(double)nblk);
         rows_in_block = rows_in_buffer_A_now;
         for(j = 0; j < intNumber; j++)
         {
            if((j+1)*nblk < cols_in_buffer_A_now)
               cols_in_block = nblk;
            else
               cols_in_block = cols_in_buffer_A_now - j*nblk;

            dlacpy("A", &rows_in_block, &cols_in_block, CopyFrom, &rows_in_block, CopyTo, &rows_in_block);

            CopyFrom = CopyFrom + rows_in_block*cols_in_block;
            CopyTo = CopyTo + nblk*(ratio*rows_in_block - nblk*(ratio-1)*ratio/2);
            rows_in_block = rows_in_block - ratio*nblk;
         }
      }
      else
      {
         if(my_prow > 0)
         {
            MPI_Sendrecv(Buf_to_send_A, (int) Size_send_A, MPI_DOUBLE, (int) where_to_send_A, (int) zero, Buf_to_receive_A, (int) Size_U_stored, MPI_DOUBLE, (int) from_where_to_receive_A, (int) zero, row_comm, &status);
            MPI_Get_count(&status, MPI_DOUBLE, &Size_receive_AMPI);
            Size_receive_A = (int) Size_receive_AMPI;

            cols_in_buffer_A = (int)Buf_to_receive_A[Size_receive_A-1];
            if(from_where_to_receive_A <= my_prow)
            {
               rows_in_buffer_A = na_rows;
            }
            else
            {
               rows_in_buffer_A = na_rows - ceil((double)(((double)from_where_to_receive_A - (double)my_prow)/(double)np_rows))*nblk;
            }
         }
         else
         {
            Size_receive_A = Size_send_A;
            rows_in_buffer_A = rows_in_buffer_A_my_initial;
            cols_in_buffer_A = cols_in_buffer_A_my_initial;
         }
      }
   }
   if(ratio != 1)
   {
      Buf_to_receive_A[Size_receive_A] = cols_in_buffer_A;
      Buf_to_receive_A[Size_receive_A + 1] = rows_in_buffer_A;
      Size_receive_A = Size_receive_A + 2;
   }
   else
   {
      Buf_to_receive_A[Size_receive_A] = rows_in_buffer_A;
      Size_receive_A = Size_receive_A + 1;
   }

   Size_receive_U = Size_U_skewed;
   U_to_calc = U_stored;

   where_to_send_A = (my_pcol - 1 + np_cols)%np_cols;
   from_where_to_receive_A = (my_pcol + 1)%np_cols;
   where_to_send_U = (my_prow - 1 + np_rows)%np_rows;
   from_where_to_receive_U = (my_prow + 1)%np_rows;
   Curr_pos_in_U_stored = Size_U_skewed;

   for(j = 1; j < np_rows; j++)
   {

      data_ptr = Buf_to_send_A;
      Buf_to_send_A = Buf_to_receive_A;
      Buf_to_receive_A = data_ptr;

      if (j > ToStore)
      {
         data_ptr = Buf_to_send_U;
         Buf_to_send_U = Buf_to_receive_U;
         Buf_to_receive_U = data_ptr;
      }

      Size_send_A = Size_receive_A;
      MPI_Isend(Buf_to_send_A, (int) Size_send_A, MPI_DOUBLE, (int) where_to_send_A, (int) zero, row_comm, &request_A_Send);
      MPI_Irecv(Buf_to_receive_A, (int) (ratio*Size_U_stored), MPI_DOUBLE, (int) from_where_to_receive_A, (int) zero, row_comm, &request_A_Recv);

      Size_send_U = Size_receive_U;
      if (j > ToStore)
      {
         if(j > ToStore + 1)
         {
            MPI_Isend(Buf_to_send_U, (int) Size_send_U, MPI_DOUBLE, (int) where_to_send_U, (int) zero, col_comm, &request_U_Send);
            U_to_calc = Buf_to_send_U;
         }
         else {
     MPI_Isend(U_to_calc, (int) Size_send_U, MPI_DOUBLE, (int) where_to_send_U, (int) zero, col_comm, &request_U_Send);
  }
         MPI_Irecv(Buf_to_receive_U, (int) Size_U_stored, MPI_DOUBLE, (int) from_where_to_receive_U, (int) zero, col_comm, &request_U_Recv);
      }

      rows_in_buffer_U = (int)U_to_calc[Size_receive_U-1];
      row_of_origin_U = (my_pcol + my_prow + np_cols + j - 1)%np_rows;
      if(my_pcol >= row_of_origin_U)
         cols_in_buffer_U = na_cols;
      else
         cols_in_buffer_U = na_cols - nblk;

      cols_in_buffer_A = (int)Buf_to_send_A[Size_receive_A-2];
      rows_in_buffer_A = (int)Buf_to_send_A[Size_receive_A-1];

      col_of_origin_A = np_cols;
      for(i = 0; i < ratio; i++)
      {
         intNumber = (my_pcol + my_prow + i*np_rows + np_cols + j - 1)%np_cols;
         if(intNumber < col_of_origin_A)
            col_of_origin_A = intNumber;
      }

      if (my_pcol >= row_of_origin_U)
         curr_col_loc_res = 0;
      else
         curr_col_loc_res = nblk;

      num_of_blocks_in_U_buffer = ceil((double)((double)cols_in_buffer_U/(double)nblk));
      if(my_pcol >= row_of_origin_U)
         rows_in_block_U = ceil(((double)(my_pcol + 1) - (double)row_of_origin_U)/(double)np_rows)*nblk;
      else
         rows_in_block_U = ratio*nblk;

      U_local_start = U_to_calc;

      for (i = 0; i < num_of_blocks_in_U_buffer; i++)
      {

         curr_col_glob_res = np_cols*nblk*(curr_col_loc_res/nblk) + curr_col_loc_res%nblk + ((np_cols+my_pcol)%np_cols)*nblk;

         Nb = curr_col_glob_res/nblk;
         owner = Nb%np_rows;
         curr_row_loc_res = (Nb/np_rows)*nblk;
         if(my_prow < owner)
            curr_row_loc_res = curr_row_loc_res + nblk;

         curr_row_loc_A = curr_row_loc_res;
         if(col_of_origin_A > my_prow)
            curr_row_loc_A = curr_row_loc_A - nblk;

         rows_in_block = rows_in_buffer_A - curr_row_loc_A;

         curr_col_loc_U = i*nblk;

         if((curr_col_loc_U + nblk) <= cols_in_buffer_U)
            cols_in_block = nblk;
         else
            cols_in_block = cols_in_buffer_U - curr_col_loc_U;

         if(rows_in_block_U > rows_in_buffer_U)
            rows_in_block_U = rows_in_buffer_U;

         A_local_index = curr_row_loc_A;
         A_local_start = &Buf_to_send_A[A_local_index];
         Res_ptr = &Res[curr_col_loc_res*na_rows + curr_row_loc_res];

         LDA_A = rows_in_buffer_A;
         LDA_A_new = LDA_A;
         if ((rows_in_block > 0)&&(cols_in_block > 0))
         {
            U_local_start_curr = U_local_start;

            for (ii = 0; ii < ceil((double)rows_in_block_U/(double)nblk); ii++)
            {
               if((ii+1)*nblk <= cols_in_buffer_A)
                  rows_in_block_U_curr = nblk;
               else
                  rows_in_block_U_curr = cols_in_buffer_A - ii*nblk;

               if((j == 1)&&(ii == 0)) {
                  dgemm_("N", "N", &rows_in_block, &cols_in_block, &rows_in_block_U_curr, &done, A_local_start, &LDA_A, U_local_start_curr, &rows_in_block_U, &dzero, Res_ptr, &na_rows);
        }
               else {
                  dgemm_("N", "N", &rows_in_block, &cols_in_block, &rows_in_block_U_curr, &done, A_local_start, &LDA_A, U_local_start_curr, &rows_in_block_U, &done, Res_ptr, &na_rows);
               }

               LDA_A_new = LDA_A_new - nblk;

               U_local_start_curr = U_local_start_curr + rows_in_block_U_curr;
               A_local_index = A_local_index - LDA_A + LDA_A*nblk + LDA_A_new;
               A_local_start = &Buf_to_send_A[A_local_index];
               LDA_A = LDA_A_new;
            }
         }

         U_local_start = U_local_start + rows_in_block_U*cols_in_block;
         curr_col_loc_res = curr_col_loc_res + nblk;
         rows_in_block_U = rows_in_block_U + ratio*nblk;
      }

      MPI_Wait(&request_A_Send, &status);
      MPI_Wait(&request_A_Recv, &status);
      MPI_Get_count(&status, MPI_DOUBLE, &Size_receive_AMPI);
      Size_receive_A = (int) Size_receive_AMPI;

      if (j <= ToStore)
      {
         U_to_calc = &U_stored[Curr_pos_in_U_stored];
         Curr_pos_in_U_stored = Curr_pos_in_U_stored + SizesU[j-1];
         Size_receive_U = SizesU[j-1];
      }
      else
      {
         MPI_Wait(&request_U_Send, &status);
         MPI_Wait(&request_U_Recv, &status);
  MPI_Get_count(&status, MPI_DOUBLE, &Size_receive_UMPI);
         Size_receive_U = (int) Size_receive_UMPI;
      }
   }

   if(ToStore < np_rows - 1)
      U_to_calc = Buf_to_receive_U;
   rows_in_buffer_U = (int)U_to_calc[Size_receive_U-1];
   row_of_origin_U = (my_pcol + my_prow + np_cols + j - 1)%np_rows;
   if(my_pcol >= row_of_origin_U)
      cols_in_buffer_U = na_cols;
   else
      cols_in_buffer_U = na_cols - nblk;

   cols_in_buffer_A = (int)Buf_to_receive_A[Size_receive_A-2];
   rows_in_buffer_A = (int)Buf_to_receive_A[Size_receive_A-1];

   col_of_origin_A = np_cols;
   for(i = 0; i < ratio; i++)
   {
      intNumber = (my_pcol + my_prow + i*np_rows + np_cols + np_rows - 1)%np_cols;
      if(intNumber < col_of_origin_A)
         col_of_origin_A = intNumber;
   }

   if (my_pcol >= row_of_origin_U)
      curr_col_loc_res = 0;
   else
      curr_col_loc_res = nblk;

   num_of_blocks_in_U_buffer = ceil((double)((double)cols_in_buffer_U/(double)nblk));
   if(my_pcol >= row_of_origin_U)
      rows_in_block_U = ceil(((double)(my_pcol + 1) - (double)row_of_origin_U)/(double)np_rows)*nblk;
   else
      rows_in_block_U = ratio*nblk;

   U_local_start = U_to_calc;

   for (i = 0; i < num_of_blocks_in_U_buffer; i++)
   {

      curr_col_glob_res = np_cols*nblk*(curr_col_loc_res/nblk) + curr_col_loc_res%nblk + ((np_cols+my_pcol)%np_cols)*nblk;

      Nb = curr_col_glob_res/nblk;
      owner = Nb%np_rows;
      curr_row_loc_res = (Nb/np_rows)*nblk;
      if(my_prow < owner)
         curr_row_loc_res = curr_row_loc_res + nblk;

      curr_row_loc_A = curr_row_loc_res;
      if(col_of_origin_A > my_prow)
         curr_row_loc_A = curr_row_loc_A - nblk;

      rows_in_block = rows_in_buffer_A - curr_row_loc_A;

      curr_col_loc_U = i*nblk;

      if((curr_col_loc_U + nblk) <= cols_in_buffer_U)
         cols_in_block = nblk;
      else
         cols_in_block = cols_in_buffer_U - curr_col_loc_U;

      if(rows_in_block_U > rows_in_buffer_U)
         rows_in_block_U = rows_in_buffer_U;

      A_local_index = curr_row_loc_A;
      A_local_start = &Buf_to_receive_A[A_local_index];
      Res_ptr = &Res[curr_col_loc_res*na_rows + curr_row_loc_res];
      LDA_A = rows_in_buffer_A;
      LDA_A_new = LDA_A;
      if ((rows_in_block > 0) &&(cols_in_block > 0))
      {
         U_local_start_curr = U_local_start;

         for (ii = 0; ii < ceil((double)rows_in_block_U/(double)nblk); ii++)
         {
            if((ii+1)*nblk <= cols_in_buffer_A)
               rows_in_block_U_curr = nblk;
            else
               rows_in_block_U_curr = cols_in_buffer_A - ii*nblk;

            if((j == 1)&&(ii == 0)) {
               dgemm_("N", "N", &rows_in_block, &cols_in_block, &rows_in_block_U_curr, &done, A_local_start, &LDA_A, U_local_start_curr, &rows_in_block_U, &dzero, Res_ptr, &na_rows);
     }
            else {
               dgemm_("N", "N", &rows_in_block, &cols_in_block, &rows_in_block_U_curr, &done, A_local_start, &LDA_A, U_local_start_curr, &rows_in_block_U, &done, Res_ptr, &na_rows);
     }

            LDA_A_new = LDA_A_new - nblk;

            U_local_start_curr = U_local_start_curr + rows_in_block_U_curr;
            A_local_index = A_local_index - (LDA_A - rows_in_block) + LDA_A*nblk + LDA_A_new - rows_in_block;
            A_local_start = &Buf_to_receive_A[A_local_index];
            LDA_A = LDA_A_new;
         }
      }

      U_local_start = U_local_start + rows_in_block_U*cols_in_block;
      curr_col_loc_res = curr_col_loc_res + nblk;
      rows_in_block_U = rows_in_block_U + ratio*nblk;
   }

   pdtran(&na, &na, &done, Res, &one, &one, a_desc, &dzero, M, &one, &one, a_desc);
   pdlacpy("U", &na, &na, M, &one, &one, a_desc, Res, &one, &one, a_desc);

   free(Buf_to_send_A);
   free(Buf_to_receive_A);
   free(Buf_to_send_U);
   free(Buf_to_receive_U);
   free(M);
   free(M_T);
   if(ratio != 1)
      free(Buf_A);
   free(U_stored);
   free(SizesU);
}

void cannons_reduction_c_d(double* A, double* U, int local_rowsCast, int local_colsCast,
                         int* a_desc, double *Res, int ToStore, int row_comm, int col_comm)
{
  int local_rows, local_cols;
  local_rows = (int) local_rowsCast;
  local_cols = (int) local_colsCast;

  MPI_Comm c_row_comm = MPI_Comm_f2c(row_comm);
  MPI_Comm c_col_comm = MPI_Comm_f2c(col_comm);

  int my_prow, my_pcol, np_rows, np_cols;
  int my_prowMPI, my_pcolMPI, np_rowsMPI, np_colsMPI;

  MPI_Comm_rank(c_row_comm, &my_prowMPI);
  MPI_Comm_size(c_row_comm, &np_rowsMPI);
  MPI_Comm_rank(c_col_comm, &my_pcolMPI);
  MPI_Comm_size(c_col_comm, &np_colsMPI);

  my_prow = (int) my_prowMPI;
  my_pcol = (int) my_pcolMPI;
  np_rows = (int) np_rowsMPI;
  np_cols = (int) np_colsMPI;

  cannons_reduction_d(A, U, np_rows, np_cols, my_prow, my_pcol, a_desc, Res, ToStore, c_col_comm, c_row_comm);
}
void cannons_triang_rectangular_d(double* U, double* B, int np_rows, int np_cols, int my_prow, int my_pcol, int* U_desc, int* b_desc, double *Res, MPI_Comm row_comm, MPI_Comm col_comm)
{
   int na, nb, nblk, width, na_rows, na_cols, nb_cols, cols_in_buffer_U_my_initial, cols_in_buffer_U, rows_in_buffer_U, Size_receive_U_now, rows_in_buffer_U_now, cols_in_buffer_U_now, rows_in_buffer_U_my_initial;

   int Size_receive_U_nowMPI, Size_receive_UMPI, Size_receive_BMPI;
   int i, j, Size_send_U, Size_receive_U, Size_send_B, Size_receive_B, intNumber, Buf_rows, Buf_cols_U, Buf_cols_B, curr_rows, num_of_iters, cols_in_buffer, rows_in_block, curr_col_loc, cols_in_block, num_of_blocks_in_U_buffer, col_of_origin_U, b_rows_mult, b_cols_mult;

   double *Buf_to_send_U, *Buf_to_receive_U, *Buf_to_send_B, *Buf_to_receive_B, *Buf_U, *PosBuff;

   int where_to_send_U, from_where_to_receive_U, where_to_send_B, from_where_to_receive_B, last_proc_col_B, last_proc_row_B, n, Size_U_stored, proc_col_min;

   double *U_local_start, *Buf_pos, *B_local_start, *double_ptr, *CopyTo, *CopyFrom;

   int ratio;

   MPI_Status status;

   int one = 1;
   int zero = 0;
   double done = 1.0;
   double dzero = 0.0;

   na = U_desc[2];
   nblk = U_desc[4];
   nb = b_desc[3];

   na_rows = numroc(&na, &nblk, &my_prow, &zero, &np_rows);
   na_cols = numroc(&na, &nblk, &my_pcol, &zero, &np_cols);
   nb_cols = numroc(&nb, &nblk, &my_pcol, &zero, &np_cols);

   MPI_Request request_U_Recv;
   MPI_Request request_U_Send;
   MPI_Request request_B_Recv;
   MPI_Request request_B_Send;

   last_proc_col_B = ((nb-1)/nblk) % np_cols;
   last_proc_row_B = ((na-1)/nblk) % np_rows;

    if(nb%nblk == 0)
      if(my_pcol <= last_proc_col_B)
         Buf_cols_B = nb_cols;
      else
         Buf_cols_B = nb_cols + nblk;
   else
      if(my_pcol < last_proc_col_B)
         Buf_cols_B = nb_cols;
      else if(my_pcol > last_proc_col_B)
         Buf_cols_B = nb_cols + nblk;
      else
         Buf_cols_B = nb_cols + nblk - nb_cols%nblk;

   if(na%nblk == 0)
      if(my_prow <= last_proc_row_B)
         Buf_rows = na_rows;
      else
         Buf_rows = na_rows + nblk;
   else
      if(my_prow < last_proc_row_B)
         Buf_rows = na_rows;
      else if(my_prow > last_proc_row_B)
         Buf_rows = na_rows + nblk;
      else
         Buf_rows = na_rows + nblk - na_rows%nblk;

   ratio = np_cols/np_rows;

   intNumber = ceil((double)na/(double)(np_cols*nblk));
   Size_U_stored = ratio*nblk*nblk*intNumber*(intNumber+1)/2 + 2;

   Buf_to_send_U = malloc(ratio*Size_U_stored*sizeof(double));
   Buf_to_receive_U = malloc(ratio*Size_U_stored*sizeof(double));
   Buf_to_send_B = malloc(Buf_cols_B*Buf_rows*sizeof(double));
   Buf_to_receive_B = malloc(Buf_cols_B*Buf_rows*sizeof(double));
   if(ratio != 1)
      Buf_U = malloc(Size_U_stored*sizeof(double));

   for(i = 0; i < na_rows*nb_cols; i++)
     Res[i] = 0;

   if((ratio != 1)||(my_prow != 0))
      Buf_pos = Buf_to_send_U;
   else
      Buf_pos = Buf_to_receive_U;

   if(my_pcol >= my_prow)
      curr_col_loc = 0;
   else
      curr_col_loc = 1;

   num_of_iters = ceil((double)na_cols/(double)nblk);
   num_of_iters = num_of_iters - curr_col_loc;
   curr_col_loc = curr_col_loc*nblk;

   if(my_pcol >= my_prow )
      rows_in_block = ceil(((double)(my_pcol + 1) - (double)my_prow)/(double)np_rows)*nblk;
   else
      rows_in_block = ratio*nblk;
   cols_in_buffer_U_my_initial = 0;
   Size_send_U = 0;
   for(i = 0; i < num_of_iters; i++)
   {
      if(rows_in_block > na_rows)
         rows_in_block = na_rows;

      if ((na_cols - curr_col_loc) < nblk)
         cols_in_block = na_cols - curr_col_loc;
      else
         cols_in_block = nblk;

      if((rows_in_block > 0)&&(cols_in_block > 0))
      {
         double_ptr = &U[curr_col_loc*na_rows];
         dlacpy("A", &rows_in_block, &cols_in_block, double_ptr, &na_rows, Buf_pos, &rows_in_block);
         Buf_pos = Buf_pos + rows_in_block*cols_in_block;
         Size_send_U = Size_send_U + rows_in_block*cols_in_block;
         cols_in_buffer_U_my_initial = cols_in_buffer_U_my_initial + cols_in_block;
      }
      curr_col_loc = curr_col_loc + nblk;
      rows_in_block = rows_in_block + ratio*nblk;
   }
   rows_in_buffer_U_my_initial = rows_in_block - ratio*nblk;
   *Buf_pos = (double)cols_in_buffer_U_my_initial;
   Buf_pos = Buf_pos + 1;
   *Buf_pos = (double)rows_in_buffer_U_my_initial;
   Size_send_U = Size_send_U + 2;

   proc_col_min = np_cols;
   for(i = 0; i < ratio; i++)
   {
      from_where_to_receive_U = (my_pcol + my_prow + i*np_rows)%np_cols;
      if(from_where_to_receive_U < proc_col_min)
         proc_col_min = from_where_to_receive_U;
   }

   Size_receive_U = 0;
   cols_in_buffer_U = 0;
   rows_in_buffer_U = 0;
   for(i = 0; i < ratio; i++)
   {
      where_to_send_U = (my_pcol - my_prow - i*np_rows + np_cols)%np_cols;
      from_where_to_receive_U = (my_pcol + my_prow + i*np_rows)%np_cols;

      if(ratio != 1)
      {
         if(where_to_send_U != my_pcol)
         {
            MPI_Sendrecv(Buf_to_send_U, (int) Size_send_U, MPI_DOUBLE, (int) where_to_send_U, 0, Buf_U, (int) Size_U_stored, MPI_DOUBLE, (int) from_where_to_receive_U, 0, row_comm, &status);
            MPI_Get_count(&status, MPI_DOUBLE, &Size_receive_U_nowMPI);
            Size_receive_U_now = (int) Size_receive_U_nowMPI;
            Size_receive_U = Size_receive_U + Size_receive_U_now - 2;

            cols_in_buffer_U_now = Buf_U[Size_receive_U_now - 2];
            cols_in_buffer_U = cols_in_buffer_U + cols_in_buffer_U_now;
            rows_in_buffer_U_now = Buf_U[Size_receive_U_now - 1];

            if(rows_in_buffer_U < rows_in_buffer_U_now)
               rows_in_buffer_U = rows_in_buffer_U_now;

            intNumber = from_where_to_receive_U/np_rows;
            if(proc_col_min >= my_prow)
               CopyTo = &Buf_to_receive_U[nblk*nblk*intNumber*(intNumber + 1)/2];
            else
               if(from_where_to_receive_U < my_prow)
                  CopyTo = &Buf_to_receive_U[nblk*nblk*ratio*(ratio - 1)/2];
               else
                  CopyTo = &Buf_to_receive_U[nblk*nblk*intNumber*(intNumber - 1)/2];
            CopyFrom = Buf_U;
         }
         else
         {
            cols_in_buffer_U_now = cols_in_buffer_U_my_initial;
            cols_in_buffer_U = cols_in_buffer_U + cols_in_buffer_U_now;

            rows_in_buffer_U_now = rows_in_buffer_U_my_initial;
            if(rows_in_buffer_U < rows_in_buffer_U_now)
               rows_in_buffer_U = rows_in_buffer_U_now;

            intNumber = my_pcol/np_rows;
            if(proc_col_min >= my_prow)
               CopyTo = &Buf_to_receive_U[nblk*nblk*intNumber*(intNumber + 1)/2];
            else
               if(my_pcol < my_prow)
                  CopyTo = &Buf_to_receive_U[nblk*nblk*ratio*(ratio - 1)/2];
               else
                  CopyTo = &Buf_to_receive_U[nblk*nblk*intNumber*(intNumber - 1)/2];
            CopyFrom = Buf_to_send_U;
            Size_receive_U = Size_receive_U + Size_send_U - 2;
         }

         intNumber = ceil((double)cols_in_buffer_U_now/(double)nblk);
         if(from_where_to_receive_U >= my_prow)
            rows_in_block = ceil(((double)(from_where_to_receive_U + 1) - (double)my_prow)/(double)np_rows)*nblk;
         else
            rows_in_block = ratio*nblk;
         for(j = 0; j < intNumber; j++)
         {
            if((j+1)*nblk < cols_in_buffer_U_now)
               cols_in_block = nblk;
            else
               cols_in_block = cols_in_buffer_U_now - j*nblk;

            dlacpy("A", &rows_in_block, &cols_in_block, CopyFrom, &rows_in_block, CopyTo, &rows_in_block);

            CopyFrom = CopyFrom + rows_in_block*cols_in_block;
            CopyTo = CopyTo + ratio*rows_in_block*nblk + nblk*nblk*ratio*(ratio-1)/2;
            rows_in_block = rows_in_block + ratio*nblk;
            if(rows_in_block > rows_in_buffer_U_now)
               rows_in_block = rows_in_buffer_U_now;
         }
      }
      else
      {
         if(my_prow > 0)
         {
            MPI_Sendrecv(Buf_to_send_U, (int) Size_send_U, MPI_DOUBLE, (int) where_to_send_U, 0, Buf_to_receive_U, (int) Size_U_stored, MPI_DOUBLE, (int) from_where_to_receive_U, 0, row_comm, &status);
            MPI_Get_count(&status, MPI_DOUBLE, &Size_receive_UMPI);
            Size_receive_U = (int) Size_receive_UMPI;

            cols_in_buffer_U = (int)Buf_to_receive_U[Size_receive_U-2];
            rows_in_buffer_U = (int)Buf_to_receive_U[Size_receive_U-1];
         }
         else
         {
            Size_receive_U = Size_send_U;
            rows_in_buffer_U = rows_in_buffer_U_my_initial;
            cols_in_buffer_U = cols_in_buffer_U_my_initial;
         }
      }
   }
   if(ratio != 1)
   {
      Buf_to_receive_U[Size_receive_U] = cols_in_buffer_U;
      Buf_to_receive_U[Size_receive_U + 1] = rows_in_buffer_U;
      Size_receive_U = Size_receive_U + 2;
   }

   if(my_pcol > 0)
   {
      where_to_send_B = (my_prow - my_pcol + np_cols)%np_rows;
      from_where_to_receive_B = (my_pcol + my_prow)%np_rows;

      if(where_to_send_B != my_prow)
      {

         dlacpy("A", &na_rows, &nb_cols, B, &na_rows, Buf_to_send_B, &na_rows);
         MPI_Sendrecv(Buf_to_send_B, (int) (nb_cols*na_rows), MPI_DOUBLE, (int) where_to_send_B, 0, Buf_to_receive_B, (int) (nb_cols*Buf_rows), MPI_DOUBLE, (int) from_where_to_receive_B, 0, col_comm, &status);
         MPI_Get_count(&status, MPI_DOUBLE, &Size_receive_BMPI);
         Size_receive_B = (int) Size_receive_BMPI;
         Size_receive_B = Size_receive_B/nb_cols;

      }
      else
      {
         dlacpy("A", &na_rows, &nb_cols, B, &na_rows, Buf_to_receive_B, &na_rows);
         Size_receive_B = na_rows;
      }
   }
   else
   {
      dlacpy("A", &na_rows, &nb_cols, B, &na_rows, Buf_to_receive_B, &na_rows);
      Size_receive_B = na_rows;
   }

   where_to_send_U = (my_pcol - 1 + np_cols)%np_cols;
   from_where_to_receive_U = (my_pcol + 1)%np_cols;
   where_to_send_B = (my_prow - 1 + np_rows)%np_rows;
   from_where_to_receive_B = (my_prow + 1)%np_rows;

   for(i = 1; i < np_rows; i++)
   {

      double_ptr = Buf_to_send_U;
      Buf_to_send_U = Buf_to_receive_U;
      Buf_to_receive_U = double_ptr;

      double_ptr = Buf_to_send_B;
      Buf_to_send_B = Buf_to_receive_B;
      Buf_to_receive_B = double_ptr;

      Size_send_U = Size_receive_U;
      Size_send_B = Size_receive_B;

      MPI_Isend(Buf_to_send_U, (int) Size_send_U, MPI_DOUBLE, (int) where_to_send_U, 0, row_comm, &request_U_Send);
      MPI_Irecv(Buf_to_receive_U, (int) (ratio*Size_U_stored), MPI_DOUBLE, (int) from_where_to_receive_U, 0, row_comm, &request_U_Recv);

      MPI_Isend(Buf_to_send_B, (int) (Size_send_B*nb_cols), MPI_DOUBLE, (int) where_to_send_B, 0, col_comm, &request_B_Send);
      MPI_Irecv(Buf_to_receive_B, (int) (Buf_rows*nb_cols), MPI_DOUBLE, (int) from_where_to_receive_B, 0, col_comm, &request_B_Recv);

      cols_in_buffer_U = (int)Buf_to_send_U[Size_receive_U-2];
      rows_in_buffer_U = (int)Buf_to_send_U[Size_receive_U-1];

      proc_col_min = np_cols;
      for(j = 0; j < ratio; j++)
      {
         col_of_origin_U = (my_pcol + my_prow + i - 1 + j*np_rows)%np_cols;
         if(col_of_origin_U < proc_col_min)
            proc_col_min = col_of_origin_U;
      }
      col_of_origin_U = proc_col_min;

      num_of_blocks_in_U_buffer = ceil((double)cols_in_buffer_U/(double)nblk);

      if (col_of_origin_U >= my_prow)
         B_local_start = Buf_to_send_B;
      else
         B_local_start = Buf_to_send_B + nblk;

      U_local_start = Buf_to_send_U;

      for(j = 0; j < num_of_blocks_in_U_buffer; j++)
      {
         curr_rows = (j+1)*nblk;
         if (curr_rows > rows_in_buffer_U)
            curr_rows = rows_in_buffer_U;

         if((j+1)*nblk <= cols_in_buffer_U)
            b_rows_mult = nblk;
         else
            b_rows_mult = cols_in_buffer_U - j*nblk;

         dgemm_("N", "N", &curr_rows, &nb_cols, &b_rows_mult, &done, U_local_start, &curr_rows, B_local_start, &Size_receive_B, &done, Res, &na_rows);

         U_local_start = U_local_start + nblk*curr_rows;
         B_local_start = B_local_start + nblk;
      }

      MPI_Wait(&request_U_Send, &status);
      MPI_Wait(&request_U_Recv, &status);
      MPI_Get_count(&status, MPI_DOUBLE, &Size_receive_UMPI);
      Size_receive_U = (int) Size_receive_UMPI;

      MPI_Wait(&request_B_Send, &status);
      MPI_Wait(&request_B_Recv, &status);
      MPI_Get_count(&status, MPI_DOUBLE, &Size_receive_BMPI);
      Size_receive_B = (int) Size_receive_BMPI;
      Size_receive_B = (int) Size_receive_B / nb_cols;

   }

   cols_in_buffer_U = (int)Buf_to_receive_U[Size_receive_U-2];
   rows_in_buffer_U = (int)Buf_to_receive_U[Size_receive_U-1];

   proc_col_min = np_cols;
   for(j = 0; j < ratio; j++)
   {
      col_of_origin_U = (my_pcol + my_prow + np_rows - 1 + j*np_rows)%np_cols;
      if(col_of_origin_U < proc_col_min)
         proc_col_min = col_of_origin_U;
   }
   col_of_origin_U = proc_col_min;

   num_of_blocks_in_U_buffer = ceil((double)cols_in_buffer_U/(double)nblk);

   if (col_of_origin_U >= my_prow)
      B_local_start = Buf_to_receive_B;
   else
      B_local_start = Buf_to_receive_B + nblk;

   U_local_start = Buf_to_receive_U;

   for(j = 0; j < num_of_blocks_in_U_buffer; j++)
   {
      curr_rows = (j+1)*nblk;
      if (curr_rows > rows_in_buffer_U)
         curr_rows = rows_in_buffer_U;

      if((j+1)*nblk <= cols_in_buffer_U)
         b_rows_mult = nblk;
      else
         b_rows_mult = cols_in_buffer_U - j*nblk;

      dgemm_("N", "N", &curr_rows, &nb_cols, &b_rows_mult, &done, U_local_start, &curr_rows, B_local_start, &Size_receive_B, &done, Res, &na_rows);

      U_local_start = U_local_start + nblk*curr_rows;
      B_local_start = B_local_start + nblk;
   }

   free(Buf_to_send_U);
   free(Buf_to_receive_U);
   free(Buf_to_send_B);
   free(Buf_to_receive_B);
   if(ratio != 1)
      free(Buf_U);
}

void cannons_triang_rectangular_c_d(double* U, double* B, int local_rowsCast, int local_colsCast,
                                    int* u_desc, int* b_desc, double *Res, int row_comm, int col_comm)
{
  int local_rows, local_cols;

  local_rows = (int) local_rowsCast;
  local_cols = (int) local_colsCast;

  MPI_Comm c_row_comm = MPI_Comm_f2c(row_comm);
  MPI_Comm c_col_comm = MPI_Comm_f2c(col_comm);

  int my_prow, my_pcol, np_rows, np_cols;
  int my_prowMPI, my_pcolMPI, np_rowsMPI, np_colsMPI;

  MPI_Comm_rank(c_row_comm, &my_prowMPI);
  MPI_Comm_size(c_row_comm, &np_rowsMPI);
  MPI_Comm_rank(c_col_comm, &my_pcolMPI);
  MPI_Comm_size(c_col_comm, &np_colsMPI);

  my_prow = (int) my_prowMPI;
  my_pcol = (int) my_pcolMPI;
  np_rows = (int) np_rowsMPI;
  np_cols = (int) np_colsMPI;

  cannons_triang_rectangular_d(U, B, np_rows, np_cols, my_prow, my_pcol, u_desc, b_desc, Res, c_col_comm, c_row_comm);
}
void cannons_reduction_c_d(double* A, double* U, int local_rowsCast, int local_colsCast, int* a_desc,
                           double *Res, int ToStore, int row_comm, int col_comm);
void cannons_triang_rectangular_c_d(double* U, double* B, int local_rowsCast, int local_colsCast,
                                    int* u_desc, int* b_desc, double *Res, int row_comm, int col_comm);

void dlacpy(char*, int*, int*, double*, int*, double*, int*);
void dgemm_(char*, char*, int*, int*, int*, double*, double*, int*, double*, int*, double*, double*, int*);

void slacpy(char*, int*, int*, float*, int*, float*, int*);
void sgemm_(char*, char*, int*, int*, int*, float*, float*, int*, float*, int*, float*, float*, int*);

void zlacpy(char*, int*, int*, double complex*, int*, double complex*, int*);
void zgemm_(char*, char*, int*, int*, int*, double complex*, double complex*, int*, double complex*, int*, double complex*, double complex*, int*);

void clacpy(char*, int*, int*, float complex*, int*, float complex*, int*);
void cgemm_(char*, char*, int*, int*, int*, float complex*, float complex*, int*, float complex*, int*, float complex*, float complex*, int*);
int numroc(int*, int*, int*, int*, int*);

void pdlacpy(char*, int*, int*, double*, int*, int*, int*, double*, int*, int*, int*);
void pdtran(int*, int*, double*, double*, int*, int*, int*, double*, double*, int*, int*, int*);

void pslacpy(char*, int*, int*, float*, int*, int*, int*, float*, int*, int*, int*);
void pstran(int*, int*, float*, float*, int*, int*, int*, float*, float*, int*, int*, int*);

void pzlacpy(char*, int*, int*, double complex*, int*, int*, int*, double complex*, int*, int*, int*);
void pztranc(int*, int*, double complex*, double complex*, int*, int*, int*, double complex*, double complex*, int*, int*, int*);

void pclacpy(char*, int*, int*, float complex*, int*, int*, int*, float complex*, int*, int*, int*);
void pctranc(int*, int*, float complex*, float complex*, int*, int*, int*, float complex*, float complex*, int*, int*, int*);

void cannons_reduction_f(float* A, float* U, int np_rows, int np_cols, int my_prow, int my_pcol,
                         int* a_desc, float *Res, int ToStore, MPI_Comm row_comm, MPI_Comm col_comm)
{
   int na, nblk, i, j, Size_send_A, Size_receive_A, Size_send_U, Size_receive_U, Buf_rows, Buf_cols, where_to_send_A, from_where_to_receive_A, where_to_send_U, from_where_to_receive_U, last_proc_row, last_proc_col, cols_in_buffer_A, rows_in_buffer_A, intNumber;
   float *Buf_to_send_A, *Buf_to_receive_A, *Buf_to_send_U, *Buf_to_receive_U, *data_ptr, *Buf_A, *Buf_pos, *U_local_start, *Res_ptr, *M, *M_T, *A_local_start, *U_local_start_curr, *U_stored, *CopyTo, *CopyFrom, *U_to_calc;
   int ratio, num_of_iters, cols_in_buffer, rows_in_block, rows_in_buffer, curr_col_loc, cols_in_block, curr_col_glob, curr_row_loc, Size_receive_A_now, Nb, owner, cols_in_buffer_A_now;
   int Size_receive_A_nowMPI, Size_receive_AMPI, Size_receive_UMPI;

   int row_of_origin_U, rows_in_block_U, num_of_blocks_in_U_buffer, k, startPos, cols_in_buffer_U, rows_in_buffer_U, col_of_origin_A, curr_row_loc_res, curr_row_loc_A, curr_col_glob_res;
   int curr_col_loc_res, curr_col_loc_buf, proc_row_curr, curr_col_loc_U, A_local_index, LDA_A, LDA_A_new, index_row_A_for_LDA, ii, rows_in_block_U_curr, width, row_origin_U, rows_in_block_A, cols_in_buffer_A_my_initial, rows_in_buffer_A_my_initial, proc_col_min;
   int *SizesU;
   int Size_U_skewed, Size_U_stored, Curr_pos_in_U_stored, rows_in_buffer_A_now;
   float done = 1.0;
   float dzero = 0.0;
   int one = 1;
   int zero = 0;
   int na_rows, na_cols;

   MPI_Status status;
   MPI_Request request_A_Recv;
   MPI_Request request_A_Send;
   MPI_Request request_U_Recv;
   MPI_Request request_U_Send;

   na = a_desc[2];
   nblk = a_desc[4];
   na_rows = numroc(&na, &nblk, &my_prow, &zero, &np_rows);
   na_cols = numroc(&na, &nblk, &my_pcol, &zero, &np_cols);

   if(ToStore > (np_rows -1))
      if((my_prow == 0)&&(my_pcol == 0))
         printf("Buffering level is larger than (np_rows-1) !!!\n");
   if((my_prow == 0)&&(my_pcol == 0))
         printf("Buffering level = %d\n", ToStore);

   if (np_cols%np_rows != 0)
   {
      if((my_prow == 0)&& (my_pcol ==0))
         printf("!!!!! np_cols must be a multiple of np_rows!!!!! I do nothing! \n");
      return;
   }
   if (np_cols < np_rows != 0)
   {
      if((my_prow == 0)&& (my_pcol ==0))
         printf("np_cols < np_rows \n");
      return;
   }

   ratio = np_cols/np_rows;
   last_proc_row = ((na-1)/nblk) % np_rows;
   last_proc_col = ((na-1)/nblk) % np_cols;

   if(na%nblk == 0)
      if(my_pcol <= last_proc_col)
         Buf_cols = na_cols;
      else
         Buf_cols = na_cols + nblk;
   else
      if(my_pcol < last_proc_col)
         Buf_cols = na_cols;
      else if(my_pcol > last_proc_col)
         Buf_cols = na_cols + nblk;
      else
         Buf_cols = na_cols + nblk - na_cols%nblk;

  if(na%nblk == 0)
      if(my_prow <= last_proc_row)
         Buf_rows = na_rows;
      else
         Buf_rows = na_rows + nblk;
   else
      if(my_prow < last_proc_row)
         Buf_rows = na_rows;
      else if(my_prow > last_proc_row)
         Buf_rows = na_rows + nblk;
      else
         Buf_rows = na_rows + nblk - na_rows%nblk;

   intNumber = ceil((float)na/(float)(np_cols*nblk));
   Size_U_stored = ratio*nblk*nblk*intNumber*(intNumber+1)/2 + 2;

   U_stored = malloc((Size_U_stored*(ToStore+1))*sizeof(float));
   SizesU = malloc(ToStore*sizeof(int));
   Buf_to_send_A = malloc(ratio*Buf_cols*Buf_rows*sizeof(float));
   Buf_to_receive_A = malloc(ratio*Buf_cols*Buf_rows*sizeof(float));
   Buf_to_send_U = malloc(Size_U_stored*sizeof(float));
   Buf_to_receive_U = malloc(Size_U_stored*sizeof(float));
   if(ratio != 1)
      Buf_A = malloc(Buf_cols*Buf_rows*sizeof(float));
   M = malloc(na_rows*na_cols*sizeof(float));
   M_T = malloc(na_rows*na_cols*sizeof(float));
   for(i = 0; i < na_rows*na_cols; i++)
      M[i] = 0;

   if(ratio != 1)
      slacpy("A", &na_rows, &na_cols, A, &na_rows, Buf_to_send_A, &na_rows);
   Size_receive_A = 0;

   for(i = 0; i < ratio; i++)
   {
      where_to_send_A = (my_pcol - my_prow - i*np_rows + np_cols)%np_cols;
      from_where_to_receive_A = (my_pcol + my_prow + i*np_rows)%np_cols;

      if(ratio != 1)
      {
         if(where_to_send_A != my_pcol)
         {
           MPI_Sendrecv(Buf_to_send_A, (int) (na_cols*na_rows), MPI_FLOAT,(int) where_to_send_A, (int) zero, Buf_A, (int) (na_rows*Buf_cols), MPI_FLOAT, (int) from_where_to_receive_A, (int) zero, row_comm, &status);
           MPI_Get_count(&status, MPI_FLOAT, &Size_receive_A_nowMPI);
           Size_receive_A_now = (int) Size_receive_A_nowMPI;
           Size_receive_A_now = Size_receive_A_now/na_rows;
         }
         else
            Size_receive_A_now = na_cols;
         Size_receive_A = Size_receive_A + Size_receive_A_now;

         intNumber = from_where_to_receive_A/np_rows;

         CopyTo = &Buf_to_receive_A[intNumber*na_rows*nblk];
         if(where_to_send_A != my_pcol)
            CopyFrom = Buf_A;
         else
            CopyFrom = A;

         intNumber = ceil((float)Size_receive_A_now/(float)nblk);
         for(j = 0; j < intNumber; j++)
         {
            width = nblk;
            if(nblk*(j+1) > Size_receive_A_now)
               width = Size_receive_A_now - nblk*j;
            slacpy("A", &na_rows, &width, CopyFrom, &na_rows, CopyTo, &na_rows);
            CopyTo = CopyTo + na_rows*nblk*ratio;
            CopyFrom = CopyFrom + na_rows*nblk;
         }
      }
      else
         if(my_prow > 0)
         {
            slacpy("A", &na_rows, &na_cols, A, &na_rows, Buf_to_send_A, &na_rows);
            MPI_Sendrecv(Buf_to_send_A, (int) (na_cols*na_rows), MPI_FLOAT, (int) where_to_send_A, (int) zero, Buf_to_receive_A, (int) (na_rows*Buf_cols), MPI_FLOAT, (int) from_where_to_receive_A, (int) zero, row_comm, &status);
            MPI_Get_count(&status, MPI_FLOAT, &Size_receive_AMPI);
            Size_receive_A = (int) Size_receive_AMPI;
     Size_receive_A = Size_receive_A/na_rows;
         }
         else
         {
            slacpy("A", &na_rows, &na_cols, A, &na_rows, Buf_to_receive_A, &na_rows);
            Size_receive_A = na_cols;
         }
   }

   num_of_iters = ceil((float)na_cols/(float)nblk);

   where_to_send_U = (my_prow - my_pcol + np_cols)%np_rows;
   from_where_to_receive_U = (my_pcol + my_prow)%np_rows;

   if(where_to_send_U == my_prow)
      Buf_pos = Buf_to_receive_U;
   else
      Buf_pos = Buf_to_send_U;

   if(my_pcol >= my_prow)
      curr_col_loc = 0;
   else
      curr_col_loc = 1;

   num_of_iters = num_of_iters - curr_col_loc;
   curr_col_loc = curr_col_loc*nblk;

   if(my_pcol >= my_prow )
      rows_in_block = ceil(((float)(my_pcol + 1) - (float)my_prow)/(float)np_rows)*nblk;
   else
      rows_in_block = ratio*nblk;

   Size_send_U = 0;
   for(i = 0; i < num_of_iters; i++)
   {
      if(rows_in_block > na_rows)
         rows_in_block = na_rows;

      if ((na_cols - curr_col_loc) < nblk)
         cols_in_block = na_cols - curr_col_loc;
      else
         cols_in_block = nblk;

      if((rows_in_block > 0)&&(cols_in_block > 0))
      {
         data_ptr = &U[curr_col_loc*na_rows];
         slacpy("A", &rows_in_block, &cols_in_block, data_ptr, &na_rows, Buf_pos, &rows_in_block);
         Buf_pos = Buf_pos + rows_in_block*cols_in_block;
         Size_send_U = Size_send_U + rows_in_block*cols_in_block;
      }
      curr_col_loc = curr_col_loc + nblk;
      rows_in_block = rows_in_block + ratio*nblk;
   }
   rows_in_buffer = rows_in_block - ratio*nblk;
   *Buf_pos = (float)rows_in_buffer;
   Size_send_U = Size_send_U + 1;

   if(where_to_send_U != my_prow)
   {

      MPI_Sendrecv(Buf_to_send_U, (int) Size_send_U, MPI_FLOAT, (int) where_to_send_U, (int) zero, Buf_to_receive_U, (int) (Buf_rows*na_cols), MPI_FLOAT, (int) from_where_to_receive_U, (int) zero, col_comm, &status);
      MPI_Get_count(&status, MPI_FLOAT, &Size_receive_UMPI);
      Size_receive_U = (int) Size_receive_UMPI;
   }
   else
      Size_receive_U = Size_send_U;

   for(i = 0; i < Size_receive_U; i++)
      U_stored[i] = Buf_to_receive_U[i];
   Size_U_skewed = Size_receive_U;
   Curr_pos_in_U_stored = Size_U_skewed;

   where_to_send_A = (my_pcol - 1 + np_cols)%np_cols;
   from_where_to_receive_A = (my_pcol + 1)%np_cols;
   where_to_send_U = (my_prow - 1 + np_rows)%np_rows;
   from_where_to_receive_U = (my_prow + 1)%np_rows;

   for(j = 1; j < np_rows; j++)
   {

      data_ptr = Buf_to_send_A;
      Buf_to_send_A = Buf_to_receive_A;
      Buf_to_receive_A = data_ptr;

      data_ptr = Buf_to_send_U;
      Buf_to_send_U = Buf_to_receive_U;
      Buf_to_receive_U = data_ptr;

      Size_send_A = Size_receive_A;
      MPI_Isend(Buf_to_send_A, (int) (Size_send_A*na_rows), MPI_FLOAT, (int) where_to_send_A, (int) zero, row_comm, &request_A_Send);
      MPI_Irecv(Buf_to_receive_A, (int) (Buf_cols*na_rows*ratio), MPI_FLOAT, (int) from_where_to_receive_A, (int) zero, row_comm, &request_A_Recv);

      Size_send_U = Size_receive_U;
      MPI_Isend(Buf_to_send_U, (int) Size_send_U, MPI_FLOAT, (int) where_to_send_U, (int) zero, col_comm, &request_U_Send);
      MPI_Irecv(Buf_to_receive_U, (int) (Buf_rows*na_cols), MPI_FLOAT, (int) from_where_to_receive_U, (int) zero, col_comm, &request_U_Recv);

      rows_in_buffer = (int)Buf_to_send_U[Size_receive_U-1];
      row_origin_U = (my_pcol + my_prow + np_cols + j - 1)%np_rows;

      if((my_pcol >= my_prow)&&(my_pcol >= row_origin_U))
      {
         cols_in_buffer = na_cols;
         curr_col_loc_res = 0;
         curr_col_loc_buf = 0;
      }
      if((my_pcol < my_prow)&&(my_pcol < row_origin_U))
      {
         cols_in_buffer = na_cols - nblk;
         curr_col_loc_res = nblk;
         curr_col_loc_buf = 0;
      }
      if((my_pcol >= my_prow)&&(my_pcol < row_origin_U))
      {
         cols_in_buffer = na_cols - nblk;
         curr_col_loc_res = nblk;
         curr_col_loc_buf = 0;
      }
      if((my_pcol < my_prow)&&(my_pcol >= row_origin_U))
      {
         cols_in_buffer = na_cols;
         curr_col_loc_res = nblk;
         curr_col_loc_buf = nblk;
      }

      num_of_blocks_in_U_buffer = ceil(((float)cols_in_buffer - (float)curr_col_loc_buf)/(float)nblk);

      startPos = (curr_col_loc_buf + nblk)*curr_col_loc_buf/2;
      U_local_start = &Buf_to_send_U[startPos];
      Res_ptr = &M[curr_col_loc_res*na_rows];

      for (i = 0; i < num_of_blocks_in_U_buffer; i++)
      {
         curr_col_glob = (curr_col_loc_res/nblk)*nblk*np_cols + my_pcol*nblk;
         proc_row_curr = (curr_col_glob/nblk)%np_rows;
         rows_in_block_A = (curr_col_glob/(nblk*np_rows))*nblk;
         if(my_prow <= proc_row_curr)
            rows_in_block_A = rows_in_block_A + nblk;

         if(rows_in_block_A > na_rows)
            rows_in_block_A = na_rows;

         if((curr_col_loc_buf + nblk) <= cols_in_buffer)
            cols_in_block = nblk;
         else
            cols_in_block = cols_in_buffer - curr_col_loc_buf;

         rows_in_block_U = (curr_col_glob/(nblk*np_rows))*nblk;
         if(proc_row_curr >= row_origin_U)
            rows_in_block_U = rows_in_block_U + nblk;

         if(rows_in_block_U > rows_in_buffer)
            rows_in_block_U = rows_in_buffer;

         if ((rows_in_block_A > 0)&&(cols_in_block > 0))
            if (j == 1) {
               sgemm_("N", "N", &rows_in_block_A, &cols_in_block, &rows_in_block_U, &done, Buf_to_send_A, &na_rows, U_local_start, &rows_in_block_U, &dzero, Res_ptr, &na_rows);
     }
            else {
               sgemm_("N", "N", &rows_in_block_A, &cols_in_block, &rows_in_block_U, &done, Buf_to_send_A, &na_rows, U_local_start, &rows_in_block_U, &done, Res_ptr, &na_rows);
     }

         U_local_start = U_local_start + rows_in_block_U*cols_in_block;
         curr_col_loc_res = curr_col_loc_res + nblk;
         Res_ptr = &M[curr_col_loc_res*na_rows];
         curr_col_loc_buf = curr_col_loc_buf + nblk;
      }

      MPI_Wait(&request_A_Send, &status);
      MPI_Wait(&request_A_Recv, &status);

      MPI_Get_count(&status, MPI_FLOAT, &Size_receive_AMPI);
      Size_receive_A = (int) Size_receive_AMPI;
      Size_receive_A = Size_receive_A / na_rows;

      MPI_Wait(&request_U_Send, &status);
      MPI_Wait(&request_U_Recv, &status);
      MPI_Get_count(&status, MPI_FLOAT, &Size_receive_UMPI);
      Size_receive_U = (int) Size_receive_UMPI;

      if(j <= ToStore)
      {
         for(k = 0; k < Size_receive_U; k++)
            U_stored[Curr_pos_in_U_stored + k] = Buf_to_receive_U[k];
         Curr_pos_in_U_stored = Curr_pos_in_U_stored + Size_receive_U;
         SizesU[j-1] = Size_receive_U;
      }
   }

   rows_in_buffer = (int)Buf_to_receive_U[Size_receive_U-1];
   row_origin_U = (my_pcol + my_prow + np_cols + np_rows -1)%np_rows;

   if((my_pcol >= my_prow)&&(my_pcol >= row_origin_U))
   {
      cols_in_buffer = na_cols;
      curr_col_loc_res = 0;
      curr_col_loc_buf = 0;
   }
   if((my_pcol < my_prow)&&(my_pcol < row_origin_U))
   {
      cols_in_buffer = na_cols - nblk;
      curr_col_loc_res = nblk;
      curr_col_loc_buf = 0;
   }
   if((my_pcol >= my_prow)&&(my_pcol < row_origin_U))
   {
      cols_in_buffer = na_cols - nblk;
      curr_col_loc_res = nblk;
      curr_col_loc_buf = 0;
   }
   if((my_pcol < my_prow)&&(my_pcol >= row_origin_U))
   {
      cols_in_buffer = na_cols;
      curr_col_loc_res = nblk;
      curr_col_loc_buf = nblk;
   }

   num_of_blocks_in_U_buffer = ceil(((float)cols_in_buffer - (float)curr_col_loc_buf)/(float)nblk);

   startPos = (curr_col_loc_buf + nblk)*curr_col_loc_buf/2;
   U_local_start = &Buf_to_receive_U[startPos];
   Res_ptr = &M[curr_col_loc_res*na_rows];

   for (i = 0; i < num_of_blocks_in_U_buffer; i++)
   {
      curr_col_glob = (curr_col_loc_res/nblk)*nblk*np_cols + my_pcol*nblk;
      proc_row_curr = (curr_col_glob/nblk)%np_rows;
      rows_in_block_A = (curr_col_glob/(nblk*np_rows))*nblk;
      if(my_prow <= proc_row_curr)
         rows_in_block_A = rows_in_block_A + nblk;

      if(rows_in_block_A > na_rows)
         rows_in_block_A = na_rows;

      if((curr_col_loc_buf + nblk) <= cols_in_buffer)
         cols_in_block = nblk;
      else
         cols_in_block = cols_in_buffer - curr_col_loc_buf;

      rows_in_block_U = (curr_col_glob/(nblk*np_rows))*nblk;
      if(proc_row_curr >= row_origin_U)
         rows_in_block_U = rows_in_block_U + nblk;

      if(rows_in_block_U > rows_in_buffer)
         rows_in_block_U = rows_in_buffer;

      if ((rows_in_block_A > 0)&&(cols_in_block > 0))
         if (j == 1) {
            sgemm_("N", "N", &rows_in_block_A, &cols_in_block, &rows_in_block_U, &done, Buf_to_receive_A, &na_rows, U_local_start, &rows_in_block_U, &dzero, Res_ptr, &na_rows);
  }
         else {
            sgemm_("N", "N", &rows_in_block_A, &cols_in_block, &rows_in_block_U, &done, Buf_to_receive_A, &na_rows, U_local_start, &rows_in_block_U, &done, Res_ptr, &na_rows);
         }

      U_local_start = U_local_start + rows_in_block_U*cols_in_block;
      curr_col_loc_res = curr_col_loc_res + nblk;
      Res_ptr = &M[curr_col_loc_res*na_rows];
      curr_col_loc_buf = curr_col_loc_buf + nblk;
   }

   pstran(&na, &na, &done, M, &one, &one, a_desc, &dzero, M_T, &one, &one, a_desc);

   if((ratio != 1)||(my_prow != 0))
      Buf_pos = Buf_to_send_A;
   else
      Buf_pos = Buf_to_receive_A;

   num_of_iters = ceil((float)na_cols/(float)nblk);

   cols_in_buffer_A_my_initial = 0;
   Size_send_A = 0;

   if(my_pcol <= my_prow)
   {
      curr_row_loc = 0;
      rows_in_buffer_A_my_initial = na_rows;
   }
   else
   {
      curr_row_loc = ceil((float)(((float)my_pcol - (float)my_prow)/(float)np_rows))*nblk;
      rows_in_buffer_A_my_initial = na_rows - curr_row_loc;
   }

   for(i = 0; i < num_of_iters; i++)
   {
      curr_col_loc = i*nblk;
      rows_in_block = na_rows - curr_row_loc;

      if ((na_cols - curr_col_loc) < nblk)
         cols_in_block = na_cols - curr_col_loc;
      else
         cols_in_block = nblk;

      if((rows_in_block > 0)&&(cols_in_block > 0))
      {
         A_local_start = &M_T[curr_col_loc*na_rows + curr_row_loc];
         slacpy("A", &rows_in_block, &cols_in_block, A_local_start, &na_rows, Buf_pos, &rows_in_block);
         Buf_pos = Buf_pos + rows_in_block*cols_in_block;
         Size_send_A = Size_send_A + rows_in_block*cols_in_block;
         cols_in_buffer_A_my_initial = cols_in_buffer_A_my_initial + cols_in_block;
      }
      curr_row_loc = curr_row_loc + ratio*nblk;
   }
   *Buf_pos = (float)cols_in_buffer_A_my_initial;
   Size_send_A = Size_send_A + 1;

   proc_col_min = np_cols;
   for(i = 0; i < ratio; i++)
   {
      from_where_to_receive_A = (my_pcol + my_prow + i*np_rows)%np_cols;
      if(from_where_to_receive_A < proc_col_min)
         proc_col_min = from_where_to_receive_A;
   }

   Size_receive_A = 0;
   cols_in_buffer_A = 0;
   rows_in_buffer_A = 0;
   for(i = 0; i < ratio; i++)
   {
      where_to_send_A = (my_pcol - my_prow - i*np_rows + np_cols)%np_cols;
      from_where_to_receive_A = (my_pcol + my_prow + i*np_rows)%np_cols;

      if(ratio != 1)
      {
         if(where_to_send_A != my_pcol)
         {
            MPI_Sendrecv(Buf_to_send_A, (int) Size_send_A, MPI_FLOAT, (int) where_to_send_A, (int) zero, Buf_A, (int) Size_U_stored, MPI_FLOAT, (int) from_where_to_receive_A, (int) zero, row_comm, &status);
            MPI_Get_count(&status, MPI_FLOAT, &Size_receive_A_nowMPI);
            Size_receive_A_now = (int) Size_receive_A_nowMPI;

            Size_receive_A = Size_receive_A + Size_receive_A_now - 1;

            cols_in_buffer_A_now = Buf_A[Size_receive_A_now-1];
            cols_in_buffer_A = cols_in_buffer_A + cols_in_buffer_A_now;

            if(from_where_to_receive_A <= my_prow)
            {
               rows_in_buffer_A_now = na_rows;
            }
            else
            {
               rows_in_buffer_A_now = na_rows - ceil((float)(((float)from_where_to_receive_A - (float)my_prow)/(float)np_rows))*nblk;
            }
            if(rows_in_buffer_A < rows_in_buffer_A_now)
               rows_in_buffer_A = rows_in_buffer_A_now;

            intNumber = from_where_to_receive_A/np_rows;
            if(proc_col_min <= my_prow)
               CopyTo = &Buf_to_receive_A[nblk*(na_rows*intNumber - nblk*(intNumber-1)*intNumber/2)];
            else
               CopyTo = &Buf_to_receive_A[nblk*(na_rows*intNumber - nblk*intNumber*(intNumber+1)/2)];
            CopyFrom = Buf_A;
         }
         else
         {
            cols_in_buffer_A_now = cols_in_buffer_A_my_initial;
            cols_in_buffer_A = cols_in_buffer_A + cols_in_buffer_A_now;

            rows_in_buffer_A_now = rows_in_buffer_A_my_initial;
            if(rows_in_buffer_A < rows_in_buffer_A_now)
               rows_in_buffer_A = rows_in_buffer_A_now;

            intNumber = my_pcol/np_rows;
            if(proc_col_min <= my_prow)
               CopyTo = &Buf_to_receive_A[nblk*(na_rows*intNumber - nblk*(intNumber-1)*intNumber/2)];
            else
               CopyTo = &Buf_to_receive_A[nblk*(na_rows*intNumber - nblk*intNumber*(intNumber+1)/2)];
            CopyFrom = Buf_to_send_A;

            Size_receive_A = Size_receive_A + Size_send_A - 1;
         }

         intNumber = ceil((float)cols_in_buffer_A_now/(float)nblk);
         rows_in_block = rows_in_buffer_A_now;
         for(j = 0; j < intNumber; j++)
         {
            if((j+1)*nblk < cols_in_buffer_A_now)
               cols_in_block = nblk;
            else
               cols_in_block = cols_in_buffer_A_now - j*nblk;

            slacpy("A", &rows_in_block, &cols_in_block, CopyFrom, &rows_in_block, CopyTo, &rows_in_block);

            CopyFrom = CopyFrom + rows_in_block*cols_in_block;
            CopyTo = CopyTo + nblk*(ratio*rows_in_block - nblk*(ratio-1)*ratio/2);
            rows_in_block = rows_in_block - ratio*nblk;
         }
      }
      else
      {
         if(my_prow > 0)
         {
            MPI_Sendrecv(Buf_to_send_A, (int) Size_send_A, MPI_FLOAT, (int) where_to_send_A, (int) zero, Buf_to_receive_A, (int) Size_U_stored, MPI_FLOAT, (int) from_where_to_receive_A, (int) zero, row_comm, &status);
            MPI_Get_count(&status, MPI_FLOAT, &Size_receive_AMPI);
            Size_receive_A = (int) Size_receive_AMPI;

            cols_in_buffer_A = (int)Buf_to_receive_A[Size_receive_A-1];
            if(from_where_to_receive_A <= my_prow)
            {
               rows_in_buffer_A = na_rows;
            }
            else
            {
               rows_in_buffer_A = na_rows - ceil((float)(((float)from_where_to_receive_A - (float)my_prow)/(float)np_rows))*nblk;
            }
         }
         else
         {
            Size_receive_A = Size_send_A;
            rows_in_buffer_A = rows_in_buffer_A_my_initial;
            cols_in_buffer_A = cols_in_buffer_A_my_initial;
         }
      }
   }
   if(ratio != 1)
   {
      Buf_to_receive_A[Size_receive_A] = cols_in_buffer_A;
      Buf_to_receive_A[Size_receive_A + 1] = rows_in_buffer_A;
      Size_receive_A = Size_receive_A + 2;
   }
   else
   {
      Buf_to_receive_A[Size_receive_A] = rows_in_buffer_A;
      Size_receive_A = Size_receive_A + 1;
   }

   Size_receive_U = Size_U_skewed;
   U_to_calc = U_stored;

   where_to_send_A = (my_pcol - 1 + np_cols)%np_cols;
   from_where_to_receive_A = (my_pcol + 1)%np_cols;
   where_to_send_U = (my_prow - 1 + np_rows)%np_rows;
   from_where_to_receive_U = (my_prow + 1)%np_rows;
   Curr_pos_in_U_stored = Size_U_skewed;

   for(j = 1; j < np_rows; j++)
   {

      data_ptr = Buf_to_send_A;
      Buf_to_send_A = Buf_to_receive_A;
      Buf_to_receive_A = data_ptr;

      if (j > ToStore)
      {
         data_ptr = Buf_to_send_U;
         Buf_to_send_U = Buf_to_receive_U;
         Buf_to_receive_U = data_ptr;
      }

      Size_send_A = Size_receive_A;
      MPI_Isend(Buf_to_send_A, (int) Size_send_A, MPI_FLOAT, (int) where_to_send_A, (int) zero, row_comm, &request_A_Send);
      MPI_Irecv(Buf_to_receive_A, (int) (ratio*Size_U_stored), MPI_FLOAT, (int) from_where_to_receive_A, (int) zero, row_comm, &request_A_Recv);

      Size_send_U = Size_receive_U;
      if (j > ToStore)
      {
         if(j > ToStore + 1)
         {
            MPI_Isend(Buf_to_send_U, (int) Size_send_U, MPI_FLOAT, (int) where_to_send_U, (int) zero, col_comm, &request_U_Send);
            U_to_calc = Buf_to_send_U;
         }
         else {
     MPI_Isend(U_to_calc, (int) Size_send_U, MPI_FLOAT, (int) where_to_send_U, (int) zero, col_comm, &request_U_Send);
  }
         MPI_Irecv(Buf_to_receive_U, (int) Size_U_stored, MPI_FLOAT, (int) from_where_to_receive_U, (int) zero, col_comm, &request_U_Recv);
      }

      rows_in_buffer_U = (int)U_to_calc[Size_receive_U-1];
      row_of_origin_U = (my_pcol + my_prow + np_cols + j - 1)%np_rows;
      if(my_pcol >= row_of_origin_U)
         cols_in_buffer_U = na_cols;
      else
         cols_in_buffer_U = na_cols - nblk;

      cols_in_buffer_A = (int)Buf_to_send_A[Size_receive_A-2];
      rows_in_buffer_A = (int)Buf_to_send_A[Size_receive_A-1];

      col_of_origin_A = np_cols;
      for(i = 0; i < ratio; i++)
      {
         intNumber = (my_pcol + my_prow + i*np_rows + np_cols + j - 1)%np_cols;
         if(intNumber < col_of_origin_A)
            col_of_origin_A = intNumber;
      }

      if (my_pcol >= row_of_origin_U)
         curr_col_loc_res = 0;
      else
         curr_col_loc_res = nblk;

      num_of_blocks_in_U_buffer = ceil((float)((float)cols_in_buffer_U/(float)nblk));
      if(my_pcol >= row_of_origin_U)
         rows_in_block_U = ceil(((float)(my_pcol + 1) - (float)row_of_origin_U)/(float)np_rows)*nblk;
      else
         rows_in_block_U = ratio*nblk;

      U_local_start = U_to_calc;

      for (i = 0; i < num_of_blocks_in_U_buffer; i++)
      {

         curr_col_glob_res = np_cols*nblk*(curr_col_loc_res/nblk) + curr_col_loc_res%nblk + ((np_cols+my_pcol)%np_cols)*nblk;

         Nb = curr_col_glob_res/nblk;
         owner = Nb%np_rows;
         curr_row_loc_res = (Nb/np_rows)*nblk;
         if(my_prow < owner)
            curr_row_loc_res = curr_row_loc_res + nblk;

         curr_row_loc_A = curr_row_loc_res;
         if(col_of_origin_A > my_prow)
            curr_row_loc_A = curr_row_loc_A - nblk;

         rows_in_block = rows_in_buffer_A - curr_row_loc_A;

         curr_col_loc_U = i*nblk;

         if((curr_col_loc_U + nblk) <= cols_in_buffer_U)
            cols_in_block = nblk;
         else
            cols_in_block = cols_in_buffer_U - curr_col_loc_U;

         if(rows_in_block_U > rows_in_buffer_U)
            rows_in_block_U = rows_in_buffer_U;

         A_local_index = curr_row_loc_A;
         A_local_start = &Buf_to_send_A[A_local_index];
         Res_ptr = &Res[curr_col_loc_res*na_rows + curr_row_loc_res];

         LDA_A = rows_in_buffer_A;
         LDA_A_new = LDA_A;
         if ((rows_in_block > 0)&&(cols_in_block > 0))
         {
            U_local_start_curr = U_local_start;

            for (ii = 0; ii < ceil((float)rows_in_block_U/(float)nblk); ii++)
            {
               if((ii+1)*nblk <= cols_in_buffer_A)
                  rows_in_block_U_curr = nblk;
               else
                  rows_in_block_U_curr = cols_in_buffer_A - ii*nblk;

               if((j == 1)&&(ii == 0)) {
                  sgemm_("N", "N", &rows_in_block, &cols_in_block, &rows_in_block_U_curr, &done, A_local_start, &LDA_A, U_local_start_curr, &rows_in_block_U, &dzero, Res_ptr, &na_rows);
        }
               else {
                  sgemm_("N", "N", &rows_in_block, &cols_in_block, &rows_in_block_U_curr, &done, A_local_start, &LDA_A, U_local_start_curr, &rows_in_block_U, &done, Res_ptr, &na_rows);
               }

               LDA_A_new = LDA_A_new - nblk;

               U_local_start_curr = U_local_start_curr + rows_in_block_U_curr;
               A_local_index = A_local_index - LDA_A + LDA_A*nblk + LDA_A_new;
               A_local_start = &Buf_to_send_A[A_local_index];
               LDA_A = LDA_A_new;
            }
         }

         U_local_start = U_local_start + rows_in_block_U*cols_in_block;
         curr_col_loc_res = curr_col_loc_res + nblk;
         rows_in_block_U = rows_in_block_U + ratio*nblk;
      }

      MPI_Wait(&request_A_Send, &status);
      MPI_Wait(&request_A_Recv, &status);
      MPI_Get_count(&status, MPI_FLOAT, &Size_receive_AMPI);
      Size_receive_A = (int) Size_receive_AMPI;

      if (j <= ToStore)
      {
         U_to_calc = &U_stored[Curr_pos_in_U_stored];
         Curr_pos_in_U_stored = Curr_pos_in_U_stored + SizesU[j-1];
         Size_receive_U = SizesU[j-1];
      }
      else
      {
         MPI_Wait(&request_U_Send, &status);
         MPI_Wait(&request_U_Recv, &status);
  MPI_Get_count(&status, MPI_FLOAT, &Size_receive_UMPI);
         Size_receive_U = (int) Size_receive_UMPI;
      }
   }

   if(ToStore < np_rows - 1)
      U_to_calc = Buf_to_receive_U;
   rows_in_buffer_U = (int)U_to_calc[Size_receive_U-1];
   row_of_origin_U = (my_pcol + my_prow + np_cols + j - 1)%np_rows;
   if(my_pcol >= row_of_origin_U)
      cols_in_buffer_U = na_cols;
   else
      cols_in_buffer_U = na_cols - nblk;

   cols_in_buffer_A = (int)Buf_to_receive_A[Size_receive_A-2];
   rows_in_buffer_A = (int)Buf_to_receive_A[Size_receive_A-1];

   col_of_origin_A = np_cols;
   for(i = 0; i < ratio; i++)
   {
      intNumber = (my_pcol + my_prow + i*np_rows + np_cols + np_rows - 1)%np_cols;
      if(intNumber < col_of_origin_A)
         col_of_origin_A = intNumber;
   }

   if (my_pcol >= row_of_origin_U)
      curr_col_loc_res = 0;
   else
      curr_col_loc_res = nblk;

   num_of_blocks_in_U_buffer = ceil((float)((float)cols_in_buffer_U/(float)nblk));
   if(my_pcol >= row_of_origin_U)
      rows_in_block_U = ceil(((float)(my_pcol + 1) - (float)row_of_origin_U)/(float)np_rows)*nblk;
   else
      rows_in_block_U = ratio*nblk;

   U_local_start = U_to_calc;

   for (i = 0; i < num_of_blocks_in_U_buffer; i++)
   {

      curr_col_glob_res = np_cols*nblk*(curr_col_loc_res/nblk) + curr_col_loc_res%nblk + ((np_cols+my_pcol)%np_cols)*nblk;

      Nb = curr_col_glob_res/nblk;
      owner = Nb%np_rows;
      curr_row_loc_res = (Nb/np_rows)*nblk;
      if(my_prow < owner)
         curr_row_loc_res = curr_row_loc_res + nblk;

      curr_row_loc_A = curr_row_loc_res;
      if(col_of_origin_A > my_prow)
         curr_row_loc_A = curr_row_loc_A - nblk;

      rows_in_block = rows_in_buffer_A - curr_row_loc_A;

      curr_col_loc_U = i*nblk;

      if((curr_col_loc_U + nblk) <= cols_in_buffer_U)
         cols_in_block = nblk;
      else
         cols_in_block = cols_in_buffer_U - curr_col_loc_U;

      if(rows_in_block_U > rows_in_buffer_U)
         rows_in_block_U = rows_in_buffer_U;

      A_local_index = curr_row_loc_A;
      A_local_start = &Buf_to_receive_A[A_local_index];
      Res_ptr = &Res[curr_col_loc_res*na_rows + curr_row_loc_res];
      LDA_A = rows_in_buffer_A;
      LDA_A_new = LDA_A;
      if ((rows_in_block > 0) &&(cols_in_block > 0))
      {
         U_local_start_curr = U_local_start;

         for (ii = 0; ii < ceil((float)rows_in_block_U/(float)nblk); ii++)
         {
            if((ii+1)*nblk <= cols_in_buffer_A)
               rows_in_block_U_curr = nblk;
            else
               rows_in_block_U_curr = cols_in_buffer_A - ii*nblk;

            if((j == 1)&&(ii == 0)) {
               sgemm_("N", "N", &rows_in_block, &cols_in_block, &rows_in_block_U_curr, &done, A_local_start, &LDA_A, U_local_start_curr, &rows_in_block_U, &dzero, Res_ptr, &na_rows);
     }
            else {
               sgemm_("N", "N", &rows_in_block, &cols_in_block, &rows_in_block_U_curr, &done, A_local_start, &LDA_A, U_local_start_curr, &rows_in_block_U, &done, Res_ptr, &na_rows);
     }

            LDA_A_new = LDA_A_new - nblk;

            U_local_start_curr = U_local_start_curr + rows_in_block_U_curr;
            A_local_index = A_local_index - (LDA_A - rows_in_block) + LDA_A*nblk + LDA_A_new - rows_in_block;
            A_local_start = &Buf_to_receive_A[A_local_index];
            LDA_A = LDA_A_new;
         }
      }

      U_local_start = U_local_start + rows_in_block_U*cols_in_block;
      curr_col_loc_res = curr_col_loc_res + nblk;
      rows_in_block_U = rows_in_block_U + ratio*nblk;
   }

   pstran(&na, &na, &done, Res, &one, &one, a_desc, &dzero, M, &one, &one, a_desc);
   pslacpy("U", &na, &na, M, &one, &one, a_desc, Res, &one, &one, a_desc);

   free(Buf_to_send_A);
   free(Buf_to_receive_A);
   free(Buf_to_send_U);
   free(Buf_to_receive_U);
   free(M);
   free(M_T);
   if(ratio != 1)
      free(Buf_A);
   free(U_stored);
   free(SizesU);
}

void cannons_reduction_c_f(float* A, float* U, int local_rowsCast, int local_colsCast,
                         int* a_desc, float *Res, int ToStore, int row_comm, int col_comm)
{
  int local_rows, local_cols;
  local_rows = (int) local_rowsCast;
  local_cols = (int) local_colsCast;

  MPI_Comm c_row_comm = MPI_Comm_f2c(row_comm);
  MPI_Comm c_col_comm = MPI_Comm_f2c(col_comm);

  int my_prow, my_pcol, np_rows, np_cols;
  int my_prowMPI, my_pcolMPI, np_rowsMPI, np_colsMPI;

  MPI_Comm_rank(c_row_comm, &my_prowMPI);
  MPI_Comm_size(c_row_comm, &np_rowsMPI);
  MPI_Comm_rank(c_col_comm, &my_pcolMPI);
  MPI_Comm_size(c_col_comm, &np_colsMPI);

  my_prow = (int) my_prowMPI;
  my_pcol = (int) my_pcolMPI;
  np_rows = (int) np_rowsMPI;
  np_cols = (int) np_colsMPI;

  cannons_reduction_f(A, U, np_rows, np_cols, my_prow, my_pcol, a_desc, Res, ToStore, c_col_comm, c_row_comm);
}
void cannons_triang_rectangular_f(float* U, float* B, int np_rows, int np_cols, int my_prow, int my_pcol, int* U_desc, int* b_desc, float *Res, MPI_Comm row_comm, MPI_Comm col_comm)
{
   int na, nb, nblk, width, na_rows, na_cols, nb_cols, cols_in_buffer_U_my_initial, cols_in_buffer_U, rows_in_buffer_U, Size_receive_U_now, rows_in_buffer_U_now, cols_in_buffer_U_now, rows_in_buffer_U_my_initial;

   int Size_receive_U_nowMPI, Size_receive_UMPI, Size_receive_BMPI;
   int i, j, Size_send_U, Size_receive_U, Size_send_B, Size_receive_B, intNumber, Buf_rows, Buf_cols_U, Buf_cols_B, curr_rows, num_of_iters, cols_in_buffer, rows_in_block, curr_col_loc, cols_in_block, num_of_blocks_in_U_buffer, col_of_origin_U, b_rows_mult, b_cols_mult;

   float *Buf_to_send_U, *Buf_to_receive_U, *Buf_to_send_B, *Buf_to_receive_B, *Buf_U, *PosBuff;

   int where_to_send_U, from_where_to_receive_U, where_to_send_B, from_where_to_receive_B, last_proc_col_B, last_proc_row_B, n, Size_U_stored, proc_col_min;

   float *U_local_start, *Buf_pos, *B_local_start, *double_ptr, *CopyTo, *CopyFrom;

   int ratio;

   MPI_Status status;

   int one = 1;
   int zero = 0;
   float done = 1.0;
   float dzero = 0.0;

   na = U_desc[2];
   nblk = U_desc[4];
   nb = b_desc[3];

   na_rows = numroc(&na, &nblk, &my_prow, &zero, &np_rows);
   na_cols = numroc(&na, &nblk, &my_pcol, &zero, &np_cols);
   nb_cols = numroc(&nb, &nblk, &my_pcol, &zero, &np_cols);

   MPI_Request request_U_Recv;
   MPI_Request request_U_Send;
   MPI_Request request_B_Recv;
   MPI_Request request_B_Send;

   last_proc_col_B = ((nb-1)/nblk) % np_cols;
   last_proc_row_B = ((na-1)/nblk) % np_rows;

    if(nb%nblk == 0)
      if(my_pcol <= last_proc_col_B)
         Buf_cols_B = nb_cols;
      else
         Buf_cols_B = nb_cols + nblk;
   else
      if(my_pcol < last_proc_col_B)
         Buf_cols_B = nb_cols;
      else if(my_pcol > last_proc_col_B)
         Buf_cols_B = nb_cols + nblk;
      else
         Buf_cols_B = nb_cols + nblk - nb_cols%nblk;

   if(na%nblk == 0)
      if(my_prow <= last_proc_row_B)
         Buf_rows = na_rows;
      else
         Buf_rows = na_rows + nblk;
   else
      if(my_prow < last_proc_row_B)
         Buf_rows = na_rows;
      else if(my_prow > last_proc_row_B)
         Buf_rows = na_rows + nblk;
      else
         Buf_rows = na_rows + nblk - na_rows%nblk;

   ratio = np_cols/np_rows;

   intNumber = ceil((float)na/(float)(np_cols*nblk));
   Size_U_stored = ratio*nblk*nblk*intNumber*(intNumber+1)/2 + 2;

   Buf_to_send_U = malloc(ratio*Size_U_stored*sizeof(float));
   Buf_to_receive_U = malloc(ratio*Size_U_stored*sizeof(float));
   Buf_to_send_B = malloc(Buf_cols_B*Buf_rows*sizeof(float));
   Buf_to_receive_B = malloc(Buf_cols_B*Buf_rows*sizeof(float));
   if(ratio != 1)
      Buf_U = malloc(Size_U_stored*sizeof(float));

   for(i = 0; i < na_rows*nb_cols; i++)
     Res[i] = 0;

   if((ratio != 1)||(my_prow != 0))
      Buf_pos = Buf_to_send_U;
   else
      Buf_pos = Buf_to_receive_U;

   if(my_pcol >= my_prow)
      curr_col_loc = 0;
   else
      curr_col_loc = 1;

   num_of_iters = ceil((float)na_cols/(float)nblk);
   num_of_iters = num_of_iters - curr_col_loc;
   curr_col_loc = curr_col_loc*nblk;

   if(my_pcol >= my_prow )
      rows_in_block = ceil(((float)(my_pcol + 1) - (float)my_prow)/(float)np_rows)*nblk;
   else
      rows_in_block = ratio*nblk;
   cols_in_buffer_U_my_initial = 0;
   Size_send_U = 0;
   for(i = 0; i < num_of_iters; i++)
   {
      if(rows_in_block > na_rows)
         rows_in_block = na_rows;

      if ((na_cols - curr_col_loc) < nblk)
         cols_in_block = na_cols - curr_col_loc;
      else
         cols_in_block = nblk;

      if((rows_in_block > 0)&&(cols_in_block > 0))
      {
         double_ptr = &U[curr_col_loc*na_rows];
         slacpy("A", &rows_in_block, &cols_in_block, double_ptr, &na_rows, Buf_pos, &rows_in_block);
         Buf_pos = Buf_pos + rows_in_block*cols_in_block;
         Size_send_U = Size_send_U + rows_in_block*cols_in_block;
         cols_in_buffer_U_my_initial = cols_in_buffer_U_my_initial + cols_in_block;
      }
      curr_col_loc = curr_col_loc + nblk;
      rows_in_block = rows_in_block + ratio*nblk;
   }
   rows_in_buffer_U_my_initial = rows_in_block - ratio*nblk;
   *Buf_pos = (float)cols_in_buffer_U_my_initial;
   Buf_pos = Buf_pos + 1;
   *Buf_pos = (float)rows_in_buffer_U_my_initial;
   Size_send_U = Size_send_U + 2;

   proc_col_min = np_cols;
   for(i = 0; i < ratio; i++)
   {
      from_where_to_receive_U = (my_pcol + my_prow + i*np_rows)%np_cols;
      if(from_where_to_receive_U < proc_col_min)
         proc_col_min = from_where_to_receive_U;
   }

   Size_receive_U = 0;
   cols_in_buffer_U = 0;
   rows_in_buffer_U = 0;
   for(i = 0; i < ratio; i++)
   {
      where_to_send_U = (my_pcol - my_prow - i*np_rows + np_cols)%np_cols;
      from_where_to_receive_U = (my_pcol + my_prow + i*np_rows)%np_cols;

      if(ratio != 1)
      {
         if(where_to_send_U != my_pcol)
         {
            MPI_Sendrecv(Buf_to_send_U, (int) Size_send_U, MPI_FLOAT, (int) where_to_send_U, 0, Buf_U, (int) Size_U_stored, MPI_FLOAT, (int) from_where_to_receive_U, 0, row_comm, &status);
            MPI_Get_count(&status, MPI_FLOAT, &Size_receive_U_nowMPI);
            Size_receive_U_now = (int) Size_receive_U_nowMPI;
            Size_receive_U = Size_receive_U + Size_receive_U_now - 2;

            cols_in_buffer_U_now = Buf_U[Size_receive_U_now - 2];
            cols_in_buffer_U = cols_in_buffer_U + cols_in_buffer_U_now;
            rows_in_buffer_U_now = Buf_U[Size_receive_U_now - 1];

            if(rows_in_buffer_U < rows_in_buffer_U_now)
               rows_in_buffer_U = rows_in_buffer_U_now;

            intNumber = from_where_to_receive_U/np_rows;
            if(proc_col_min >= my_prow)
               CopyTo = &Buf_to_receive_U[nblk*nblk*intNumber*(intNumber + 1)/2];
            else
               if(from_where_to_receive_U < my_prow)
                  CopyTo = &Buf_to_receive_U[nblk*nblk*ratio*(ratio - 1)/2];
               else
                  CopyTo = &Buf_to_receive_U[nblk*nblk*intNumber*(intNumber - 1)/2];
            CopyFrom = Buf_U;
         }
         else
         {
            cols_in_buffer_U_now = cols_in_buffer_U_my_initial;
            cols_in_buffer_U = cols_in_buffer_U + cols_in_buffer_U_now;

            rows_in_buffer_U_now = rows_in_buffer_U_my_initial;
            if(rows_in_buffer_U < rows_in_buffer_U_now)
               rows_in_buffer_U = rows_in_buffer_U_now;

            intNumber = my_pcol/np_rows;
            if(proc_col_min >= my_prow)
               CopyTo = &Buf_to_receive_U[nblk*nblk*intNumber*(intNumber + 1)/2];
            else
               if(my_pcol < my_prow)
                  CopyTo = &Buf_to_receive_U[nblk*nblk*ratio*(ratio - 1)/2];
               else
                  CopyTo = &Buf_to_receive_U[nblk*nblk*intNumber*(intNumber - 1)/2];
            CopyFrom = Buf_to_send_U;
            Size_receive_U = Size_receive_U + Size_send_U - 2;
         }

         intNumber = ceil((float)cols_in_buffer_U_now/(float)nblk);
         if(from_where_to_receive_U >= my_prow)
            rows_in_block = ceil(((float)(from_where_to_receive_U + 1) - (float)my_prow)/(float)np_rows)*nblk;
         else
            rows_in_block = ratio*nblk;
         for(j = 0; j < intNumber; j++)
         {
            if((j+1)*nblk < cols_in_buffer_U_now)
               cols_in_block = nblk;
            else
               cols_in_block = cols_in_buffer_U_now - j*nblk;

            slacpy("A", &rows_in_block, &cols_in_block, CopyFrom, &rows_in_block, CopyTo, &rows_in_block);

            CopyFrom = CopyFrom + rows_in_block*cols_in_block;
            CopyTo = CopyTo + ratio*rows_in_block*nblk + nblk*nblk*ratio*(ratio-1)/2;
            rows_in_block = rows_in_block + ratio*nblk;
            if(rows_in_block > rows_in_buffer_U_now)
               rows_in_block = rows_in_buffer_U_now;
         }
      }
      else
      {
         if(my_prow > 0)
         {
            MPI_Sendrecv(Buf_to_send_U, (int) Size_send_U, MPI_FLOAT, (int) where_to_send_U, 0, Buf_to_receive_U, (int) Size_U_stored, MPI_FLOAT, (int) from_where_to_receive_U, 0, row_comm, &status);
            MPI_Get_count(&status, MPI_FLOAT, &Size_receive_UMPI);
            Size_receive_U = (int) Size_receive_UMPI;

            cols_in_buffer_U = (int)Buf_to_receive_U[Size_receive_U-2];
            rows_in_buffer_U = (int)Buf_to_receive_U[Size_receive_U-1];
         }
         else
         {
            Size_receive_U = Size_send_U;
            rows_in_buffer_U = rows_in_buffer_U_my_initial;
            cols_in_buffer_U = cols_in_buffer_U_my_initial;
         }
      }
   }
   if(ratio != 1)
   {
      Buf_to_receive_U[Size_receive_U] = cols_in_buffer_U;
      Buf_to_receive_U[Size_receive_U + 1] = rows_in_buffer_U;
      Size_receive_U = Size_receive_U + 2;
   }

   if(my_pcol > 0)
   {
      where_to_send_B = (my_prow - my_pcol + np_cols)%np_rows;
      from_where_to_receive_B = (my_pcol + my_prow)%np_rows;

      if(where_to_send_B != my_prow)
      {

         slacpy("A", &na_rows, &nb_cols, B, &na_rows, Buf_to_send_B, &na_rows);
         MPI_Sendrecv(Buf_to_send_B, (int) (nb_cols*na_rows), MPI_FLOAT, (int) where_to_send_B, 0, Buf_to_receive_B, (int) (nb_cols*Buf_rows), MPI_FLOAT, (int) from_where_to_receive_B, 0, col_comm, &status);
         MPI_Get_count(&status, MPI_FLOAT, &Size_receive_BMPI);
         Size_receive_B = (int) Size_receive_BMPI;
         Size_receive_B = Size_receive_B/nb_cols;

      }
      else
      {
         slacpy("A", &na_rows, &nb_cols, B, &na_rows, Buf_to_receive_B, &na_rows);
         Size_receive_B = na_rows;
      }
   }
   else
   {
      slacpy("A", &na_rows, &nb_cols, B, &na_rows, Buf_to_receive_B, &na_rows);
      Size_receive_B = na_rows;
   }

   where_to_send_U = (my_pcol - 1 + np_cols)%np_cols;
   from_where_to_receive_U = (my_pcol + 1)%np_cols;
   where_to_send_B = (my_prow - 1 + np_rows)%np_rows;
   from_where_to_receive_B = (my_prow + 1)%np_rows;

   for(i = 1; i < np_rows; i++)
   {

      double_ptr = Buf_to_send_U;
      Buf_to_send_U = Buf_to_receive_U;
      Buf_to_receive_U = double_ptr;

      double_ptr = Buf_to_send_B;
      Buf_to_send_B = Buf_to_receive_B;
      Buf_to_receive_B = double_ptr;

      Size_send_U = Size_receive_U;
      Size_send_B = Size_receive_B;

      MPI_Isend(Buf_to_send_U, (int) Size_send_U, MPI_FLOAT, (int) where_to_send_U, 0, row_comm, &request_U_Send);
      MPI_Irecv(Buf_to_receive_U, (int) (ratio*Size_U_stored), MPI_FLOAT, (int) from_where_to_receive_U, 0, row_comm, &request_U_Recv);

      MPI_Isend(Buf_to_send_B, (int) (Size_send_B*nb_cols), MPI_FLOAT, (int) where_to_send_B, 0, col_comm, &request_B_Send);
      MPI_Irecv(Buf_to_receive_B, (int) (Buf_rows*nb_cols), MPI_FLOAT, (int) from_where_to_receive_B, 0, col_comm, &request_B_Recv);

      cols_in_buffer_U = (int)Buf_to_send_U[Size_receive_U-2];
      rows_in_buffer_U = (int)Buf_to_send_U[Size_receive_U-1];

      proc_col_min = np_cols;
      for(j = 0; j < ratio; j++)
      {
         col_of_origin_U = (my_pcol + my_prow + i - 1 + j*np_rows)%np_cols;
         if(col_of_origin_U < proc_col_min)
            proc_col_min = col_of_origin_U;
      }
      col_of_origin_U = proc_col_min;

      num_of_blocks_in_U_buffer = ceil((float)cols_in_buffer_U/(float)nblk);

      if (col_of_origin_U >= my_prow)
         B_local_start = Buf_to_send_B;
      else
         B_local_start = Buf_to_send_B + nblk;

      U_local_start = Buf_to_send_U;

      for(j = 0; j < num_of_blocks_in_U_buffer; j++)
      {
         curr_rows = (j+1)*nblk;
         if (curr_rows > rows_in_buffer_U)
            curr_rows = rows_in_buffer_U;

         if((j+1)*nblk <= cols_in_buffer_U)
            b_rows_mult = nblk;
         else
            b_rows_mult = cols_in_buffer_U - j*nblk;

         sgemm_("N", "N", &curr_rows, &nb_cols, &b_rows_mult, &done, U_local_start, &curr_rows, B_local_start, &Size_receive_B, &done, Res, &na_rows);

         U_local_start = U_local_start + nblk*curr_rows;
         B_local_start = B_local_start + nblk;
      }

      MPI_Wait(&request_U_Send, &status);
      MPI_Wait(&request_U_Recv, &status);
      MPI_Get_count(&status, MPI_FLOAT, &Size_receive_UMPI);
      Size_receive_U = (int) Size_receive_UMPI;

      MPI_Wait(&request_B_Send, &status);
      MPI_Wait(&request_B_Recv, &status);
      MPI_Get_count(&status, MPI_FLOAT, &Size_receive_BMPI);
      Size_receive_B = (int) Size_receive_BMPI;
      Size_receive_B = (int) Size_receive_B / nb_cols;

   }

   cols_in_buffer_U = (int)Buf_to_receive_U[Size_receive_U-2];
   rows_in_buffer_U = (int)Buf_to_receive_U[Size_receive_U-1];

   proc_col_min = np_cols;
   for(j = 0; j < ratio; j++)
   {
      col_of_origin_U = (my_pcol + my_prow + np_rows - 1 + j*np_rows)%np_cols;
      if(col_of_origin_U < proc_col_min)
         proc_col_min = col_of_origin_U;
   }
   col_of_origin_U = proc_col_min;

   num_of_blocks_in_U_buffer = ceil((float)cols_in_buffer_U/(float)nblk);

   if (col_of_origin_U >= my_prow)
      B_local_start = Buf_to_receive_B;
   else
      B_local_start = Buf_to_receive_B + nblk;

   U_local_start = Buf_to_receive_U;

   for(j = 0; j < num_of_blocks_in_U_buffer; j++)
   {
      curr_rows = (j+1)*nblk;
      if (curr_rows > rows_in_buffer_U)
         curr_rows = rows_in_buffer_U;

      if((j+1)*nblk <= cols_in_buffer_U)
         b_rows_mult = nblk;
      else
         b_rows_mult = cols_in_buffer_U - j*nblk;

      sgemm_("N", "N", &curr_rows, &nb_cols, &b_rows_mult, &done, U_local_start, &curr_rows, B_local_start, &Size_receive_B, &done, Res, &na_rows);

      U_local_start = U_local_start + nblk*curr_rows;
      B_local_start = B_local_start + nblk;
   }

   free(Buf_to_send_U);
   free(Buf_to_receive_U);
   free(Buf_to_send_B);
   free(Buf_to_receive_B);
   if(ratio != 1)
      free(Buf_U);
}

void cannons_triang_rectangular_c_f(float* U, float* B, int local_rowsCast, int local_colsCast,
                                    int* u_desc, int* b_desc, float *Res, int row_comm, int col_comm)
{
  int local_rows, local_cols;

  local_rows = (int) local_rowsCast;
  local_cols = (int) local_colsCast;

  MPI_Comm c_row_comm = MPI_Comm_f2c(row_comm);
  MPI_Comm c_col_comm = MPI_Comm_f2c(col_comm);

  int my_prow, my_pcol, np_rows, np_cols;
  int my_prowMPI, my_pcolMPI, np_rowsMPI, np_colsMPI;

  MPI_Comm_rank(c_row_comm, &my_prowMPI);
  MPI_Comm_size(c_row_comm, &np_rowsMPI);
  MPI_Comm_rank(c_col_comm, &my_pcolMPI);
  MPI_Comm_size(c_col_comm, &np_colsMPI);

  my_prow = (int) my_prowMPI;
  my_pcol = (int) my_pcolMPI;
  np_rows = (int) np_rowsMPI;
  np_cols = (int) np_colsMPI;

  cannons_triang_rectangular_f(U, B, np_rows, np_cols, my_prow, my_pcol, u_desc, b_desc, Res, c_col_comm, c_row_comm);
}
void cannons_reduction_c_f(float* A, float* U, int local_rowsCast, int local_colsCast, int* a_desc,
                           float *Res, int ToStore, int row_comm, int col_comm);
void cannons_triang_rectangular_c_f(float* U, float* B, int local_rowsCast, int local_colsCast,
                                    int* u_desc, int* b_desc, float *Res, int row_comm, int col_comm);

void dlacpy(char*, int*, int*, double*, int*, double*, int*);
void dgemm_(char*, char*, int*, int*, int*, double*, double*, int*, double*, int*, double*, double*, int*);

void slacpy(char*, int*, int*, float*, int*, float*, int*);
void sgemm_(char*, char*, int*, int*, int*, float*, float*, int*, float*, int*, float*, float*, int*);

void zlacpy(char*, int*, int*, double complex*, int*, double complex*, int*);
void zgemm_(char*, char*, int*, int*, int*, double complex*, double complex*, int*, double complex*, int*, double complex*, double complex*, int*);

void clacpy(char*, int*, int*, float complex*, int*, float complex*, int*);
void cgemm_(char*, char*, int*, int*, int*, float complex*, float complex*, int*, float complex*, int*, float complex*, float complex*, int*);
int numroc(int*, int*, int*, int*, int*);

void pdlacpy(char*, int*, int*, double*, int*, int*, int*, double*, int*, int*, int*);
void pdtran(int*, int*, double*, double*, int*, int*, int*, double*, double*, int*, int*, int*);

void pslacpy(char*, int*, int*, float*, int*, int*, int*, float*, int*, int*, int*);
void pstran(int*, int*, float*, float*, int*, int*, int*, float*, float*, int*, int*, int*);

void pzlacpy(char*, int*, int*, double complex*, int*, int*, int*, double complex*, int*, int*, int*);
void pztranc(int*, int*, double complex*, double complex*, int*, int*, int*, double complex*, double complex*, int*, int*, int*);

void pclacpy(char*, int*, int*, float complex*, int*, int*, int*, float complex*, int*, int*, int*);
void pctranc(int*, int*, float complex*, float complex*, int*, int*, int*, float complex*, float complex*, int*, int*, int*);

void cannons_reduction_dc(double complex* A, double complex* U, int np_rows, int np_cols, int my_prow, int my_pcol,
                         int* a_desc, double complex *Res, int ToStore, MPI_Comm row_comm, MPI_Comm col_comm)
{
   int na, nblk, i, j, Size_send_A, Size_receive_A, Size_send_U, Size_receive_U, Buf_rows, Buf_cols, where_to_send_A, from_where_to_receive_A, where_to_send_U, from_where_to_receive_U, last_proc_row, last_proc_col, cols_in_buffer_A, rows_in_buffer_A, intNumber;
   double complex *Buf_to_send_A, *Buf_to_receive_A, *Buf_to_send_U, *Buf_to_receive_U, *data_ptr, *Buf_A, *Buf_pos, *U_local_start, *Res_ptr, *M, *M_T, *A_local_start, *U_local_start_curr, *U_stored, *CopyTo, *CopyFrom, *U_to_calc;
   int ratio, num_of_iters, cols_in_buffer, rows_in_block, rows_in_buffer, curr_col_loc, cols_in_block, curr_col_glob, curr_row_loc, Size_receive_A_now, Nb, owner, cols_in_buffer_A_now;
   int Size_receive_A_nowMPI, Size_receive_AMPI, Size_receive_UMPI;

   int row_of_origin_U, rows_in_block_U, num_of_blocks_in_U_buffer, k, startPos, cols_in_buffer_U, rows_in_buffer_U, col_of_origin_A, curr_row_loc_res, curr_row_loc_A, curr_col_glob_res;
   int curr_col_loc_res, curr_col_loc_buf, proc_row_curr, curr_col_loc_U, A_local_index, LDA_A, LDA_A_new, index_row_A_for_LDA, ii, rows_in_block_U_curr, width, row_origin_U, rows_in_block_A, cols_in_buffer_A_my_initial, rows_in_buffer_A_my_initial, proc_col_min;
   int *SizesU;
   int Size_U_skewed, Size_U_stored, Curr_pos_in_U_stored, rows_in_buffer_A_now;
   double complex done = 1.0;
   double complex dzero = 0.0;
   int one = 1;
   int zero = 0;
   int na_rows, na_cols;

   MPI_Status status;
   MPI_Request request_A_Recv;
   MPI_Request request_A_Send;
   MPI_Request request_U_Recv;
   MPI_Request request_U_Send;

   na = a_desc[2];
   nblk = a_desc[4];
   na_rows = numroc(&na, &nblk, &my_prow, &zero, &np_rows);
   na_cols = numroc(&na, &nblk, &my_pcol, &zero, &np_cols);

   if(ToStore > (np_rows -1))
      if((my_prow == 0)&&(my_pcol == 0))
         printf("Buffering level is larger than (np_rows-1) !!!\n");
   if((my_prow == 0)&&(my_pcol == 0))
         printf("Buffering level = %d\n", ToStore);

   if (np_cols%np_rows != 0)
   {
      if((my_prow == 0)&& (my_pcol ==0))
         printf("!!!!! np_cols must be a multiple of np_rows!!!!! I do nothing! \n");
      return;
   }
   if (np_cols < np_rows != 0)
   {
      if((my_prow == 0)&& (my_pcol ==0))
         printf("np_cols < np_rows \n");
      return;
   }

   ratio = np_cols/np_rows;
   last_proc_row = ((na-1)/nblk) % np_rows;
   last_proc_col = ((na-1)/nblk) % np_cols;

   if(na%nblk == 0)
      if(my_pcol <= last_proc_col)
         Buf_cols = na_cols;
      else
         Buf_cols = na_cols + nblk;
   else
      if(my_pcol < last_proc_col)
         Buf_cols = na_cols;
      else if(my_pcol > last_proc_col)
         Buf_cols = na_cols + nblk;
      else
         Buf_cols = na_cols + nblk - na_cols%nblk;

  if(na%nblk == 0)
      if(my_prow <= last_proc_row)
         Buf_rows = na_rows;
      else
         Buf_rows = na_rows + nblk;
   else
      if(my_prow < last_proc_row)
         Buf_rows = na_rows;
      else if(my_prow > last_proc_row)
         Buf_rows = na_rows + nblk;
      else
         Buf_rows = na_rows + nblk - na_rows%nblk;

   intNumber = ceil((double complex)na/(double complex)(np_cols*nblk));
   Size_U_stored = ratio*nblk*nblk*intNumber*(intNumber+1)/2 + 2;

   U_stored = malloc((Size_U_stored*(ToStore+1))*sizeof(double complex));
   SizesU = malloc(ToStore*sizeof(int));
   Buf_to_send_A = malloc(ratio*Buf_cols*Buf_rows*sizeof(double complex));
   Buf_to_receive_A = malloc(ratio*Buf_cols*Buf_rows*sizeof(double complex));
   Buf_to_send_U = malloc(Size_U_stored*sizeof(double complex));
   Buf_to_receive_U = malloc(Size_U_stored*sizeof(double complex));
   if(ratio != 1)
      Buf_A = malloc(Buf_cols*Buf_rows*sizeof(double complex));
   M = malloc(na_rows*na_cols*sizeof(double complex));
   M_T = malloc(na_rows*na_cols*sizeof(double complex));
   for(i = 0; i < na_rows*na_cols; i++)
      M[i] = 0;

   if(ratio != 1)
      zlacpy("A", &na_rows, &na_cols, A, &na_rows, Buf_to_send_A, &na_rows);
   Size_receive_A = 0;

   for(i = 0; i < ratio; i++)
   {
      where_to_send_A = (my_pcol - my_prow - i*np_rows + np_cols)%np_cols;
      from_where_to_receive_A = (my_pcol + my_prow + i*np_rows)%np_cols;

      if(ratio != 1)
      {
         if(where_to_send_A != my_pcol)
         {
           MPI_Sendrecv(Buf_to_send_A, (int) (na_cols*na_rows), MPI_DOUBLE_COMPLEX,(int) where_to_send_A, (int) zero, Buf_A, (int) (na_rows*Buf_cols), MPI_DOUBLE_COMPLEX, (int) from_where_to_receive_A, (int) zero, row_comm, &status);
           MPI_Get_count(&status, MPI_DOUBLE_COMPLEX, &Size_receive_A_nowMPI);
           Size_receive_A_now = (int) Size_receive_A_nowMPI;
           Size_receive_A_now = Size_receive_A_now/na_rows;
         }
         else
            Size_receive_A_now = na_cols;
         Size_receive_A = Size_receive_A + Size_receive_A_now;

         intNumber = from_where_to_receive_A/np_rows;

         CopyTo = &Buf_to_receive_A[intNumber*na_rows*nblk];
         if(where_to_send_A != my_pcol)
            CopyFrom = Buf_A;
         else
            CopyFrom = A;

         intNumber = ceil((double complex)Size_receive_A_now/(double complex)nblk);
         for(j = 0; j < intNumber; j++)
         {
            width = nblk;
            if(nblk*(j+1) > Size_receive_A_now)
               width = Size_receive_A_now - nblk*j;
            zlacpy("A", &na_rows, &width, CopyFrom, &na_rows, CopyTo, &na_rows);
            CopyTo = CopyTo + na_rows*nblk*ratio;
            CopyFrom = CopyFrom + na_rows*nblk;
         }
      }
      else
         if(my_prow > 0)
         {
            zlacpy("A", &na_rows, &na_cols, A, &na_rows, Buf_to_send_A, &na_rows);
            MPI_Sendrecv(Buf_to_send_A, (int) (na_cols*na_rows), MPI_DOUBLE_COMPLEX, (int) where_to_send_A, (int) zero, Buf_to_receive_A, (int) (na_rows*Buf_cols), MPI_DOUBLE_COMPLEX, (int) from_where_to_receive_A, (int) zero, row_comm, &status);
            MPI_Get_count(&status, MPI_DOUBLE_COMPLEX, &Size_receive_AMPI);
            Size_receive_A = (int) Size_receive_AMPI;
     Size_receive_A = Size_receive_A/na_rows;
         }
         else
         {
            zlacpy("A", &na_rows, &na_cols, A, &na_rows, Buf_to_receive_A, &na_rows);
            Size_receive_A = na_cols;
         }
   }

   num_of_iters = ceil((double complex)na_cols/(double complex)nblk);

   where_to_send_U = (my_prow - my_pcol + np_cols)%np_rows;
   from_where_to_receive_U = (my_pcol + my_prow)%np_rows;

   if(where_to_send_U == my_prow)
      Buf_pos = Buf_to_receive_U;
   else
      Buf_pos = Buf_to_send_U;

   if(my_pcol >= my_prow)
      curr_col_loc = 0;
   else
      curr_col_loc = 1;

   num_of_iters = num_of_iters - curr_col_loc;
   curr_col_loc = curr_col_loc*nblk;

   if(my_pcol >= my_prow )
      rows_in_block = ceil(((double complex)(my_pcol + 1) - (double complex)my_prow)/(double complex)np_rows)*nblk;
   else
      rows_in_block = ratio*nblk;

   Size_send_U = 0;
   for(i = 0; i < num_of_iters; i++)
   {
      if(rows_in_block > na_rows)
         rows_in_block = na_rows;

      if ((na_cols - curr_col_loc) < nblk)
         cols_in_block = na_cols - curr_col_loc;
      else
         cols_in_block = nblk;

      if((rows_in_block > 0)&&(cols_in_block > 0))
      {
         data_ptr = &U[curr_col_loc*na_rows];
         zlacpy("A", &rows_in_block, &cols_in_block, data_ptr, &na_rows, Buf_pos, &rows_in_block);
         Buf_pos = Buf_pos + rows_in_block*cols_in_block;
         Size_send_U = Size_send_U + rows_in_block*cols_in_block;
      }
      curr_col_loc = curr_col_loc + nblk;
      rows_in_block = rows_in_block + ratio*nblk;
   }
   rows_in_buffer = rows_in_block - ratio*nblk;
   *Buf_pos = (double complex)rows_in_buffer;
   Size_send_U = Size_send_U + 1;

   if(where_to_send_U != my_prow)
   {

      MPI_Sendrecv(Buf_to_send_U, (int) Size_send_U, MPI_DOUBLE_COMPLEX, (int) where_to_send_U, (int) zero, Buf_to_receive_U, (int) (Buf_rows*na_cols), MPI_DOUBLE_COMPLEX, (int) from_where_to_receive_U, (int) zero, col_comm, &status);
      MPI_Get_count(&status, MPI_DOUBLE_COMPLEX, &Size_receive_UMPI);
      Size_receive_U = (int) Size_receive_UMPI;
   }
   else
      Size_receive_U = Size_send_U;

   for(i = 0; i < Size_receive_U; i++)
      U_stored[i] = Buf_to_receive_U[i];
   Size_U_skewed = Size_receive_U;
   Curr_pos_in_U_stored = Size_U_skewed;

   where_to_send_A = (my_pcol - 1 + np_cols)%np_cols;
   from_where_to_receive_A = (my_pcol + 1)%np_cols;
   where_to_send_U = (my_prow - 1 + np_rows)%np_rows;
   from_where_to_receive_U = (my_prow + 1)%np_rows;

   for(j = 1; j < np_rows; j++)
   {

      data_ptr = Buf_to_send_A;
      Buf_to_send_A = Buf_to_receive_A;
      Buf_to_receive_A = data_ptr;

      data_ptr = Buf_to_send_U;
      Buf_to_send_U = Buf_to_receive_U;
      Buf_to_receive_U = data_ptr;

      Size_send_A = Size_receive_A;
      MPI_Isend(Buf_to_send_A, (int) (Size_send_A*na_rows), MPI_DOUBLE_COMPLEX, (int) where_to_send_A, (int) zero, row_comm, &request_A_Send);
      MPI_Irecv(Buf_to_receive_A, (int) (Buf_cols*na_rows*ratio), MPI_DOUBLE_COMPLEX, (int) from_where_to_receive_A, (int) zero, row_comm, &request_A_Recv);

      Size_send_U = Size_receive_U;
      MPI_Isend(Buf_to_send_U, (int) Size_send_U, MPI_DOUBLE_COMPLEX, (int) where_to_send_U, (int) zero, col_comm, &request_U_Send);
      MPI_Irecv(Buf_to_receive_U, (int) (Buf_rows*na_cols), MPI_DOUBLE_COMPLEX, (int) from_where_to_receive_U, (int) zero, col_comm, &request_U_Recv);

      rows_in_buffer = (int)Buf_to_send_U[Size_receive_U-1];
      row_origin_U = (my_pcol + my_prow + np_cols + j - 1)%np_rows;

      if((my_pcol >= my_prow)&&(my_pcol >= row_origin_U))
      {
         cols_in_buffer = na_cols;
         curr_col_loc_res = 0;
         curr_col_loc_buf = 0;
      }
      if((my_pcol < my_prow)&&(my_pcol < row_origin_U))
      {
         cols_in_buffer = na_cols - nblk;
         curr_col_loc_res = nblk;
         curr_col_loc_buf = 0;
      }
      if((my_pcol >= my_prow)&&(my_pcol < row_origin_U))
      {
         cols_in_buffer = na_cols - nblk;
         curr_col_loc_res = nblk;
         curr_col_loc_buf = 0;
      }
      if((my_pcol < my_prow)&&(my_pcol >= row_origin_U))
      {
         cols_in_buffer = na_cols;
         curr_col_loc_res = nblk;
         curr_col_loc_buf = nblk;
      }

      num_of_blocks_in_U_buffer = ceil(((double complex)cols_in_buffer - (double complex)curr_col_loc_buf)/(double complex)nblk);

      startPos = (curr_col_loc_buf + nblk)*curr_col_loc_buf/2;
      U_local_start = &Buf_to_send_U[startPos];
      Res_ptr = &M[curr_col_loc_res*na_rows];

      for (i = 0; i < num_of_blocks_in_U_buffer; i++)
      {
         curr_col_glob = (curr_col_loc_res/nblk)*nblk*np_cols + my_pcol*nblk;
         proc_row_curr = (curr_col_glob/nblk)%np_rows;
         rows_in_block_A = (curr_col_glob/(nblk*np_rows))*nblk;
         if(my_prow <= proc_row_curr)
            rows_in_block_A = rows_in_block_A + nblk;

         if(rows_in_block_A > na_rows)
            rows_in_block_A = na_rows;

         if((curr_col_loc_buf + nblk) <= cols_in_buffer)
            cols_in_block = nblk;
         else
            cols_in_block = cols_in_buffer - curr_col_loc_buf;

         rows_in_block_U = (curr_col_glob/(nblk*np_rows))*nblk;
         if(proc_row_curr >= row_origin_U)
            rows_in_block_U = rows_in_block_U + nblk;

         if(rows_in_block_U > rows_in_buffer)
            rows_in_block_U = rows_in_buffer;

         if ((rows_in_block_A > 0)&&(cols_in_block > 0))
            if (j == 1) {
               zgemm_("N", "N", &rows_in_block_A, &cols_in_block, &rows_in_block_U, &done, Buf_to_send_A, &na_rows, U_local_start, &rows_in_block_U, &dzero, Res_ptr, &na_rows);
     }
            else {
               zgemm_("N", "N", &rows_in_block_A, &cols_in_block, &rows_in_block_U, &done, Buf_to_send_A, &na_rows, U_local_start, &rows_in_block_U, &done, Res_ptr, &na_rows);
     }

         U_local_start = U_local_start + rows_in_block_U*cols_in_block;
         curr_col_loc_res = curr_col_loc_res + nblk;
         Res_ptr = &M[curr_col_loc_res*na_rows];
         curr_col_loc_buf = curr_col_loc_buf + nblk;
      }

      MPI_Wait(&request_A_Send, &status);
      MPI_Wait(&request_A_Recv, &status);

      MPI_Get_count(&status, MPI_DOUBLE_COMPLEX, &Size_receive_AMPI);
      Size_receive_A = (int) Size_receive_AMPI;
      Size_receive_A = Size_receive_A / na_rows;

      MPI_Wait(&request_U_Send, &status);
      MPI_Wait(&request_U_Recv, &status);
      MPI_Get_count(&status, MPI_DOUBLE_COMPLEX, &Size_receive_UMPI);
      Size_receive_U = (int) Size_receive_UMPI;

      if(j <= ToStore)
      {
         for(k = 0; k < Size_receive_U; k++)
            U_stored[Curr_pos_in_U_stored + k] = Buf_to_receive_U[k];
         Curr_pos_in_U_stored = Curr_pos_in_U_stored + Size_receive_U;
         SizesU[j-1] = Size_receive_U;
      }
   }

   rows_in_buffer = (int)Buf_to_receive_U[Size_receive_U-1];
   row_origin_U = (my_pcol + my_prow + np_cols + np_rows -1)%np_rows;

   if((my_pcol >= my_prow)&&(my_pcol >= row_origin_U))
   {
      cols_in_buffer = na_cols;
      curr_col_loc_res = 0;
      curr_col_loc_buf = 0;
   }
   if((my_pcol < my_prow)&&(my_pcol < row_origin_U))
   {
      cols_in_buffer = na_cols - nblk;
      curr_col_loc_res = nblk;
      curr_col_loc_buf = 0;
   }
   if((my_pcol >= my_prow)&&(my_pcol < row_origin_U))
   {
      cols_in_buffer = na_cols - nblk;
      curr_col_loc_res = nblk;
      curr_col_loc_buf = 0;
   }
   if((my_pcol < my_prow)&&(my_pcol >= row_origin_U))
   {
      cols_in_buffer = na_cols;
      curr_col_loc_res = nblk;
      curr_col_loc_buf = nblk;
   }

   num_of_blocks_in_U_buffer = ceil(((double complex)cols_in_buffer - (double complex)curr_col_loc_buf)/(double complex)nblk);

   startPos = (curr_col_loc_buf + nblk)*curr_col_loc_buf/2;
   U_local_start = &Buf_to_receive_U[startPos];
   Res_ptr = &M[curr_col_loc_res*na_rows];

   for (i = 0; i < num_of_blocks_in_U_buffer; i++)
   {
      curr_col_glob = (curr_col_loc_res/nblk)*nblk*np_cols + my_pcol*nblk;
      proc_row_curr = (curr_col_glob/nblk)%np_rows;
      rows_in_block_A = (curr_col_glob/(nblk*np_rows))*nblk;
      if(my_prow <= proc_row_curr)
         rows_in_block_A = rows_in_block_A + nblk;

      if(rows_in_block_A > na_rows)
         rows_in_block_A = na_rows;

      if((curr_col_loc_buf + nblk) <= cols_in_buffer)
         cols_in_block = nblk;
      else
         cols_in_block = cols_in_buffer - curr_col_loc_buf;

      rows_in_block_U = (curr_col_glob/(nblk*np_rows))*nblk;
      if(proc_row_curr >= row_origin_U)
         rows_in_block_U = rows_in_block_U + nblk;

      if(rows_in_block_U > rows_in_buffer)
         rows_in_block_U = rows_in_buffer;

      if ((rows_in_block_A > 0)&&(cols_in_block > 0))
         if (j == 1) {
            zgemm_("N", "N", &rows_in_block_A, &cols_in_block, &rows_in_block_U, &done, Buf_to_receive_A, &na_rows, U_local_start, &rows_in_block_U, &dzero, Res_ptr, &na_rows);
  }
         else {
            zgemm_("N", "N", &rows_in_block_A, &cols_in_block, &rows_in_block_U, &done, Buf_to_receive_A, &na_rows, U_local_start, &rows_in_block_U, &done, Res_ptr, &na_rows);
         }

      U_local_start = U_local_start + rows_in_block_U*cols_in_block;
      curr_col_loc_res = curr_col_loc_res + nblk;
      Res_ptr = &M[curr_col_loc_res*na_rows];
      curr_col_loc_buf = curr_col_loc_buf + nblk;
   }

   pztranc(&na, &na, &done, M, &one, &one, a_desc, &dzero, M_T, &one, &one, a_desc);

   if((ratio != 1)||(my_prow != 0))
      Buf_pos = Buf_to_send_A;
   else
      Buf_pos = Buf_to_receive_A;

   num_of_iters = ceil((double complex)na_cols/(double complex)nblk);

   cols_in_buffer_A_my_initial = 0;
   Size_send_A = 0;

   if(my_pcol <= my_prow)
   {
      curr_row_loc = 0;
      rows_in_buffer_A_my_initial = na_rows;
   }
   else
   {
      curr_row_loc = ceil((double complex)(((double complex)my_pcol - (double complex)my_prow)/(double complex)np_rows))*nblk;
      rows_in_buffer_A_my_initial = na_rows - curr_row_loc;
   }

   for(i = 0; i < num_of_iters; i++)
   {
      curr_col_loc = i*nblk;
      rows_in_block = na_rows - curr_row_loc;

      if ((na_cols - curr_col_loc) < nblk)
         cols_in_block = na_cols - curr_col_loc;
      else
         cols_in_block = nblk;

      if((rows_in_block > 0)&&(cols_in_block > 0))
      {
         A_local_start = &M_T[curr_col_loc*na_rows + curr_row_loc];
         zlacpy("A", &rows_in_block, &cols_in_block, A_local_start, &na_rows, Buf_pos, &rows_in_block);
         Buf_pos = Buf_pos + rows_in_block*cols_in_block;
         Size_send_A = Size_send_A + rows_in_block*cols_in_block;
         cols_in_buffer_A_my_initial = cols_in_buffer_A_my_initial + cols_in_block;
      }
      curr_row_loc = curr_row_loc + ratio*nblk;
   }
   *Buf_pos = (double complex)cols_in_buffer_A_my_initial;
   Size_send_A = Size_send_A + 1;

   proc_col_min = np_cols;
   for(i = 0; i < ratio; i++)
   {
      from_where_to_receive_A = (my_pcol + my_prow + i*np_rows)%np_cols;
      if(from_where_to_receive_A < proc_col_min)
         proc_col_min = from_where_to_receive_A;
   }

   Size_receive_A = 0;
   cols_in_buffer_A = 0;
   rows_in_buffer_A = 0;
   for(i = 0; i < ratio; i++)
   {
      where_to_send_A = (my_pcol - my_prow - i*np_rows + np_cols)%np_cols;
      from_where_to_receive_A = (my_pcol + my_prow + i*np_rows)%np_cols;

      if(ratio != 1)
      {
         if(where_to_send_A != my_pcol)
         {
            MPI_Sendrecv(Buf_to_send_A, (int) Size_send_A, MPI_DOUBLE_COMPLEX, (int) where_to_send_A, (int) zero, Buf_A, (int) Size_U_stored, MPI_DOUBLE_COMPLEX, (int) from_where_to_receive_A, (int) zero, row_comm, &status);
            MPI_Get_count(&status, MPI_DOUBLE_COMPLEX, &Size_receive_A_nowMPI);
            Size_receive_A_now = (int) Size_receive_A_nowMPI;

            Size_receive_A = Size_receive_A + Size_receive_A_now - 1;

            cols_in_buffer_A_now = Buf_A[Size_receive_A_now-1];
            cols_in_buffer_A = cols_in_buffer_A + cols_in_buffer_A_now;

            if(from_where_to_receive_A <= my_prow)
            {
               rows_in_buffer_A_now = na_rows;
            }
            else
            {
               rows_in_buffer_A_now = na_rows - ceil((double complex)(((double complex)from_where_to_receive_A - (double complex)my_prow)/(double complex)np_rows))*nblk;
            }
            if(rows_in_buffer_A < rows_in_buffer_A_now)
               rows_in_buffer_A = rows_in_buffer_A_now;

            intNumber = from_where_to_receive_A/np_rows;
            if(proc_col_min <= my_prow)
               CopyTo = &Buf_to_receive_A[nblk*(na_rows*intNumber - nblk*(intNumber-1)*intNumber/2)];
            else
               CopyTo = &Buf_to_receive_A[nblk*(na_rows*intNumber - nblk*intNumber*(intNumber+1)/2)];
            CopyFrom = Buf_A;
         }
         else
         {
            cols_in_buffer_A_now = cols_in_buffer_A_my_initial;
            cols_in_buffer_A = cols_in_buffer_A + cols_in_buffer_A_now;

            rows_in_buffer_A_now = rows_in_buffer_A_my_initial;
            if(rows_in_buffer_A < rows_in_buffer_A_now)
               rows_in_buffer_A = rows_in_buffer_A_now;

            intNumber = my_pcol/np_rows;
            if(proc_col_min <= my_prow)
               CopyTo = &Buf_to_receive_A[nblk*(na_rows*intNumber - nblk*(intNumber-1)*intNumber/2)];
            else
               CopyTo = &Buf_to_receive_A[nblk*(na_rows*intNumber - nblk*intNumber*(intNumber+1)/2)];
            CopyFrom = Buf_to_send_A;

            Size_receive_A = Size_receive_A + Size_send_A - 1;
         }

         intNumber = ceil((double complex)cols_in_buffer_A_now/(double complex)nblk);
         rows_in_block = rows_in_buffer_A_now;
         for(j = 0; j < intNumber; j++)
         {
            if((j+1)*nblk < cols_in_buffer_A_now)
               cols_in_block = nblk;
            else
               cols_in_block = cols_in_buffer_A_now - j*nblk;

            zlacpy("A", &rows_in_block, &cols_in_block, CopyFrom, &rows_in_block, CopyTo, &rows_in_block);

            CopyFrom = CopyFrom + rows_in_block*cols_in_block;
            CopyTo = CopyTo + nblk*(ratio*rows_in_block - nblk*(ratio-1)*ratio/2);
            rows_in_block = rows_in_block - ratio*nblk;
         }
      }
      else
      {
         if(my_prow > 0)
         {
            MPI_Sendrecv(Buf_to_send_A, (int) Size_send_A, MPI_DOUBLE_COMPLEX, (int) where_to_send_A, (int) zero, Buf_to_receive_A, (int) Size_U_stored, MPI_DOUBLE_COMPLEX, (int) from_where_to_receive_A, (int) zero, row_comm, &status);
            MPI_Get_count(&status, MPI_DOUBLE_COMPLEX, &Size_receive_AMPI);
            Size_receive_A = (int) Size_receive_AMPI;

            cols_in_buffer_A = (int)Buf_to_receive_A[Size_receive_A-1];
            if(from_where_to_receive_A <= my_prow)
            {
               rows_in_buffer_A = na_rows;
            }
            else
            {
               rows_in_buffer_A = na_rows - ceil((double complex)(((double complex)from_where_to_receive_A - (double complex)my_prow)/(double complex)np_rows))*nblk;
            }
         }
         else
         {
            Size_receive_A = Size_send_A;
            rows_in_buffer_A = rows_in_buffer_A_my_initial;
            cols_in_buffer_A = cols_in_buffer_A_my_initial;
         }
      }
   }
   if(ratio != 1)
   {
      Buf_to_receive_A[Size_receive_A] = cols_in_buffer_A;
      Buf_to_receive_A[Size_receive_A + 1] = rows_in_buffer_A;
      Size_receive_A = Size_receive_A + 2;
   }
   else
   {
      Buf_to_receive_A[Size_receive_A] = rows_in_buffer_A;
      Size_receive_A = Size_receive_A + 1;
   }

   Size_receive_U = Size_U_skewed;
   U_to_calc = U_stored;

   where_to_send_A = (my_pcol - 1 + np_cols)%np_cols;
   from_where_to_receive_A = (my_pcol + 1)%np_cols;
   where_to_send_U = (my_prow - 1 + np_rows)%np_rows;
   from_where_to_receive_U = (my_prow + 1)%np_rows;
   Curr_pos_in_U_stored = Size_U_skewed;

   for(j = 1; j < np_rows; j++)
   {

      data_ptr = Buf_to_send_A;
      Buf_to_send_A = Buf_to_receive_A;
      Buf_to_receive_A = data_ptr;

      if (j > ToStore)
      {
         data_ptr = Buf_to_send_U;
         Buf_to_send_U = Buf_to_receive_U;
         Buf_to_receive_U = data_ptr;
      }

      Size_send_A = Size_receive_A;
      MPI_Isend(Buf_to_send_A, (int) Size_send_A, MPI_DOUBLE_COMPLEX, (int) where_to_send_A, (int) zero, row_comm, &request_A_Send);
      MPI_Irecv(Buf_to_receive_A, (int) (ratio*Size_U_stored), MPI_DOUBLE_COMPLEX, (int) from_where_to_receive_A, (int) zero, row_comm, &request_A_Recv);

      Size_send_U = Size_receive_U;
      if (j > ToStore)
      {
         if(j > ToStore + 1)
         {
            MPI_Isend(Buf_to_send_U, (int) Size_send_U, MPI_DOUBLE_COMPLEX, (int) where_to_send_U, (int) zero, col_comm, &request_U_Send);
            U_to_calc = Buf_to_send_U;
         }
         else {
     MPI_Isend(U_to_calc, (int) Size_send_U, MPI_DOUBLE_COMPLEX, (int) where_to_send_U, (int) zero, col_comm, &request_U_Send);
  }
         MPI_Irecv(Buf_to_receive_U, (int) Size_U_stored, MPI_DOUBLE_COMPLEX, (int) from_where_to_receive_U, (int) zero, col_comm, &request_U_Recv);
      }

      rows_in_buffer_U = (int)U_to_calc[Size_receive_U-1];
      row_of_origin_U = (my_pcol + my_prow + np_cols + j - 1)%np_rows;
      if(my_pcol >= row_of_origin_U)
         cols_in_buffer_U = na_cols;
      else
         cols_in_buffer_U = na_cols - nblk;

      cols_in_buffer_A = (int)Buf_to_send_A[Size_receive_A-2];
      rows_in_buffer_A = (int)Buf_to_send_A[Size_receive_A-1];

      col_of_origin_A = np_cols;
      for(i = 0; i < ratio; i++)
      {
         intNumber = (my_pcol + my_prow + i*np_rows + np_cols + j - 1)%np_cols;
         if(intNumber < col_of_origin_A)
            col_of_origin_A = intNumber;
      }

      if (my_pcol >= row_of_origin_U)
         curr_col_loc_res = 0;
      else
         curr_col_loc_res = nblk;

      num_of_blocks_in_U_buffer = ceil((double complex)((double complex)cols_in_buffer_U/(double complex)nblk));
      if(my_pcol >= row_of_origin_U)
         rows_in_block_U = ceil(((double complex)(my_pcol + 1) - (double complex)row_of_origin_U)/(double complex)np_rows)*nblk;
      else
         rows_in_block_U = ratio*nblk;

      U_local_start = U_to_calc;

      for (i = 0; i < num_of_blocks_in_U_buffer; i++)
      {

         curr_col_glob_res = np_cols*nblk*(curr_col_loc_res/nblk) + curr_col_loc_res%nblk + ((np_cols+my_pcol)%np_cols)*nblk;

         Nb = curr_col_glob_res/nblk;
         owner = Nb%np_rows;
         curr_row_loc_res = (Nb/np_rows)*nblk;
         if(my_prow < owner)
            curr_row_loc_res = curr_row_loc_res + nblk;

         curr_row_loc_A = curr_row_loc_res;
         if(col_of_origin_A > my_prow)
            curr_row_loc_A = curr_row_loc_A - nblk;

         rows_in_block = rows_in_buffer_A - curr_row_loc_A;

         curr_col_loc_U = i*nblk;

         if((curr_col_loc_U + nblk) <= cols_in_buffer_U)
            cols_in_block = nblk;
         else
            cols_in_block = cols_in_buffer_U - curr_col_loc_U;

         if(rows_in_block_U > rows_in_buffer_U)
            rows_in_block_U = rows_in_buffer_U;

         A_local_index = curr_row_loc_A;
         A_local_start = &Buf_to_send_A[A_local_index];
         Res_ptr = &Res[curr_col_loc_res*na_rows + curr_row_loc_res];

         LDA_A = rows_in_buffer_A;
         LDA_A_new = LDA_A;
         if ((rows_in_block > 0)&&(cols_in_block > 0))
         {
            U_local_start_curr = U_local_start;

            for (ii = 0; ii < ceil((double complex)rows_in_block_U/(double complex)nblk); ii++)
            {
               if((ii+1)*nblk <= cols_in_buffer_A)
                  rows_in_block_U_curr = nblk;
               else
                  rows_in_block_U_curr = cols_in_buffer_A - ii*nblk;

               if((j == 1)&&(ii == 0)) {
                  zgemm_("N", "N", &rows_in_block, &cols_in_block, &rows_in_block_U_curr, &done, A_local_start, &LDA_A, U_local_start_curr, &rows_in_block_U, &dzero, Res_ptr, &na_rows);
        }
               else {
                  zgemm_("N", "N", &rows_in_block, &cols_in_block, &rows_in_block_U_curr, &done, A_local_start, &LDA_A, U_local_start_curr, &rows_in_block_U, &done, Res_ptr, &na_rows);
               }

               LDA_A_new = LDA_A_new - nblk;

               U_local_start_curr = U_local_start_curr + rows_in_block_U_curr;
               A_local_index = A_local_index - LDA_A + LDA_A*nblk + LDA_A_new;
               A_local_start = &Buf_to_send_A[A_local_index];
               LDA_A = LDA_A_new;
            }
         }

         U_local_start = U_local_start + rows_in_block_U*cols_in_block;
         curr_col_loc_res = curr_col_loc_res + nblk;
         rows_in_block_U = rows_in_block_U + ratio*nblk;
      }

      MPI_Wait(&request_A_Send, &status);
      MPI_Wait(&request_A_Recv, &status);
      MPI_Get_count(&status, MPI_DOUBLE_COMPLEX, &Size_receive_AMPI);
      Size_receive_A = (int) Size_receive_AMPI;

      if (j <= ToStore)
      {
         U_to_calc = &U_stored[Curr_pos_in_U_stored];
         Curr_pos_in_U_stored = Curr_pos_in_U_stored + SizesU[j-1];
         Size_receive_U = SizesU[j-1];
      }
      else
      {
         MPI_Wait(&request_U_Send, &status);
         MPI_Wait(&request_U_Recv, &status);
  MPI_Get_count(&status, MPI_DOUBLE_COMPLEX, &Size_receive_UMPI);
         Size_receive_U = (int) Size_receive_UMPI;
      }
   }

   if(ToStore < np_rows - 1)
      U_to_calc = Buf_to_receive_U;
   rows_in_buffer_U = (int)U_to_calc[Size_receive_U-1];
   row_of_origin_U = (my_pcol + my_prow + np_cols + j - 1)%np_rows;
   if(my_pcol >= row_of_origin_U)
      cols_in_buffer_U = na_cols;
   else
      cols_in_buffer_U = na_cols - nblk;

   cols_in_buffer_A = (int)Buf_to_receive_A[Size_receive_A-2];
   rows_in_buffer_A = (int)Buf_to_receive_A[Size_receive_A-1];

   col_of_origin_A = np_cols;
   for(i = 0; i < ratio; i++)
   {
      intNumber = (my_pcol + my_prow + i*np_rows + np_cols + np_rows - 1)%np_cols;
      if(intNumber < col_of_origin_A)
         col_of_origin_A = intNumber;
   }

   if (my_pcol >= row_of_origin_U)
      curr_col_loc_res = 0;
   else
      curr_col_loc_res = nblk;

   num_of_blocks_in_U_buffer = ceil((double complex)((double complex)cols_in_buffer_U/(double complex)nblk));
   if(my_pcol >= row_of_origin_U)
      rows_in_block_U = ceil(((double complex)(my_pcol + 1) - (double complex)row_of_origin_U)/(double complex)np_rows)*nblk;
   else
      rows_in_block_U = ratio*nblk;

   U_local_start = U_to_calc;

   for (i = 0; i < num_of_blocks_in_U_buffer; i++)
   {

      curr_col_glob_res = np_cols*nblk*(curr_col_loc_res/nblk) + curr_col_loc_res%nblk + ((np_cols+my_pcol)%np_cols)*nblk;

      Nb = curr_col_glob_res/nblk;
      owner = Nb%np_rows;
      curr_row_loc_res = (Nb/np_rows)*nblk;
      if(my_prow < owner)
         curr_row_loc_res = curr_row_loc_res + nblk;

      curr_row_loc_A = curr_row_loc_res;
      if(col_of_origin_A > my_prow)
         curr_row_loc_A = curr_row_loc_A - nblk;

      rows_in_block = rows_in_buffer_A - curr_row_loc_A;

      curr_col_loc_U = i*nblk;

      if((curr_col_loc_U + nblk) <= cols_in_buffer_U)
         cols_in_block = nblk;
      else
         cols_in_block = cols_in_buffer_U - curr_col_loc_U;

      if(rows_in_block_U > rows_in_buffer_U)
         rows_in_block_U = rows_in_buffer_U;

      A_local_index = curr_row_loc_A;
      A_local_start = &Buf_to_receive_A[A_local_index];
      Res_ptr = &Res[curr_col_loc_res*na_rows + curr_row_loc_res];
      LDA_A = rows_in_buffer_A;
      LDA_A_new = LDA_A;
      if ((rows_in_block > 0) &&(cols_in_block > 0))
      {
         U_local_start_curr = U_local_start;

         for (ii = 0; ii < ceil((double complex)rows_in_block_U/(double complex)nblk); ii++)
         {
            if((ii+1)*nblk <= cols_in_buffer_A)
               rows_in_block_U_curr = nblk;
            else
               rows_in_block_U_curr = cols_in_buffer_A - ii*nblk;

            if((j == 1)&&(ii == 0)) {
               zgemm_("N", "N", &rows_in_block, &cols_in_block, &rows_in_block_U_curr, &done, A_local_start, &LDA_A, U_local_start_curr, &rows_in_block_U, &dzero, Res_ptr, &na_rows);
     }
            else {
               zgemm_("N", "N", &rows_in_block, &cols_in_block, &rows_in_block_U_curr, &done, A_local_start, &LDA_A, U_local_start_curr, &rows_in_block_U, &done, Res_ptr, &na_rows);
     }

            LDA_A_new = LDA_A_new - nblk;

            U_local_start_curr = U_local_start_curr + rows_in_block_U_curr;
            A_local_index = A_local_index - (LDA_A - rows_in_block) + LDA_A*nblk + LDA_A_new - rows_in_block;
            A_local_start = &Buf_to_receive_A[A_local_index];
            LDA_A = LDA_A_new;
         }
      }

      U_local_start = U_local_start + rows_in_block_U*cols_in_block;
      curr_col_loc_res = curr_col_loc_res + nblk;
      rows_in_block_U = rows_in_block_U + ratio*nblk;
   }

   pztranc(&na, &na, &done, Res, &one, &one, a_desc, &dzero, M, &one, &one, a_desc);
   pzlacpy("U", &na, &na, M, &one, &one, a_desc, Res, &one, &one, a_desc);

   free(Buf_to_send_A);
   free(Buf_to_receive_A);
   free(Buf_to_send_U);
   free(Buf_to_receive_U);
   free(M);
   free(M_T);
   if(ratio != 1)
      free(Buf_A);
   free(U_stored);
   free(SizesU);
}

void cannons_reduction_c_dc(double complex* A, double complex* U, int local_rowsCast, int local_colsCast,
                         int* a_desc, double complex *Res, int ToStore, int row_comm, int col_comm)
{
  int local_rows, local_cols;
  local_rows = (int) local_rowsCast;
  local_cols = (int) local_colsCast;

  MPI_Comm c_row_comm = MPI_Comm_f2c(row_comm);
  MPI_Comm c_col_comm = MPI_Comm_f2c(col_comm);

  int my_prow, my_pcol, np_rows, np_cols;
  int my_prowMPI, my_pcolMPI, np_rowsMPI, np_colsMPI;

  MPI_Comm_rank(c_row_comm, &my_prowMPI);
  MPI_Comm_size(c_row_comm, &np_rowsMPI);
  MPI_Comm_rank(c_col_comm, &my_pcolMPI);
  MPI_Comm_size(c_col_comm, &np_colsMPI);

  my_prow = (int) my_prowMPI;
  my_pcol = (int) my_pcolMPI;
  np_rows = (int) np_rowsMPI;
  np_cols = (int) np_colsMPI;

  cannons_reduction_dc(A, U, np_rows, np_cols, my_prow, my_pcol, a_desc, Res, ToStore, c_col_comm, c_row_comm);
}
void cannons_triang_rectangular_dc(double complex* U, double complex* B, int np_rows, int np_cols, int my_prow, int my_pcol, int* U_desc, int* b_desc, double complex *Res, MPI_Comm row_comm, MPI_Comm col_comm)
{
   int na, nb, nblk, width, na_rows, na_cols, nb_cols, cols_in_buffer_U_my_initial, cols_in_buffer_U, rows_in_buffer_U, Size_receive_U_now, rows_in_buffer_U_now, cols_in_buffer_U_now, rows_in_buffer_U_my_initial;

   int Size_receive_U_nowMPI, Size_receive_UMPI, Size_receive_BMPI;
   int i, j, Size_send_U, Size_receive_U, Size_send_B, Size_receive_B, intNumber, Buf_rows, Buf_cols_U, Buf_cols_B, curr_rows, num_of_iters, cols_in_buffer, rows_in_block, curr_col_loc, cols_in_block, num_of_blocks_in_U_buffer, col_of_origin_U, b_rows_mult, b_cols_mult;

   double complex *Buf_to_send_U, *Buf_to_receive_U, *Buf_to_send_B, *Buf_to_receive_B, *Buf_U, *PosBuff;

   int where_to_send_U, from_where_to_receive_U, where_to_send_B, from_where_to_receive_B, last_proc_col_B, last_proc_row_B, n, Size_U_stored, proc_col_min;

   double complex *U_local_start, *Buf_pos, *B_local_start, *double_ptr, *CopyTo, *CopyFrom;

   int ratio;

   MPI_Status status;

   int one = 1;
   int zero = 0;
   double complex done = 1.0;
   double complex dzero = 0.0;

   na = U_desc[2];
   nblk = U_desc[4];
   nb = b_desc[3];

   na_rows = numroc(&na, &nblk, &my_prow, &zero, &np_rows);
   na_cols = numroc(&na, &nblk, &my_pcol, &zero, &np_cols);
   nb_cols = numroc(&nb, &nblk, &my_pcol, &zero, &np_cols);

   MPI_Request request_U_Recv;
   MPI_Request request_U_Send;
   MPI_Request request_B_Recv;
   MPI_Request request_B_Send;

   last_proc_col_B = ((nb-1)/nblk) % np_cols;
   last_proc_row_B = ((na-1)/nblk) % np_rows;

    if(nb%nblk == 0)
      if(my_pcol <= last_proc_col_B)
         Buf_cols_B = nb_cols;
      else
         Buf_cols_B = nb_cols + nblk;
   else
      if(my_pcol < last_proc_col_B)
         Buf_cols_B = nb_cols;
      else if(my_pcol > last_proc_col_B)
         Buf_cols_B = nb_cols + nblk;
      else
         Buf_cols_B = nb_cols + nblk - nb_cols%nblk;

   if(na%nblk == 0)
      if(my_prow <= last_proc_row_B)
         Buf_rows = na_rows;
      else
         Buf_rows = na_rows + nblk;
   else
      if(my_prow < last_proc_row_B)
         Buf_rows = na_rows;
      else if(my_prow > last_proc_row_B)
         Buf_rows = na_rows + nblk;
      else
         Buf_rows = na_rows + nblk - na_rows%nblk;

   ratio = np_cols/np_rows;

   intNumber = ceil((double complex)na/(double complex)(np_cols*nblk));
   Size_U_stored = ratio*nblk*nblk*intNumber*(intNumber+1)/2 + 2;

   Buf_to_send_U = malloc(ratio*Size_U_stored*sizeof(double complex));
   Buf_to_receive_U = malloc(ratio*Size_U_stored*sizeof(double complex));
   Buf_to_send_B = malloc(Buf_cols_B*Buf_rows*sizeof(double complex));
   Buf_to_receive_B = malloc(Buf_cols_B*Buf_rows*sizeof(double complex));
   if(ratio != 1)
      Buf_U = malloc(Size_U_stored*sizeof(double complex));

   for(i = 0; i < na_rows*nb_cols; i++)
     Res[i] = 0;

   if((ratio != 1)||(my_prow != 0))
      Buf_pos = Buf_to_send_U;
   else
      Buf_pos = Buf_to_receive_U;

   if(my_pcol >= my_prow)
      curr_col_loc = 0;
   else
      curr_col_loc = 1;

   num_of_iters = ceil((double complex)na_cols/(double complex)nblk);
   num_of_iters = num_of_iters - curr_col_loc;
   curr_col_loc = curr_col_loc*nblk;

   if(my_pcol >= my_prow )
      rows_in_block = ceil(((double complex)(my_pcol + 1) - (double complex)my_prow)/(double complex)np_rows)*nblk;
   else
      rows_in_block = ratio*nblk;
   cols_in_buffer_U_my_initial = 0;
   Size_send_U = 0;
   for(i = 0; i < num_of_iters; i++)
   {
      if(rows_in_block > na_rows)
         rows_in_block = na_rows;

      if ((na_cols - curr_col_loc) < nblk)
         cols_in_block = na_cols - curr_col_loc;
      else
         cols_in_block = nblk;

      if((rows_in_block > 0)&&(cols_in_block > 0))
      {
         double_ptr = &U[curr_col_loc*na_rows];
         zlacpy("A", &rows_in_block, &cols_in_block, double_ptr, &na_rows, Buf_pos, &rows_in_block);
         Buf_pos = Buf_pos + rows_in_block*cols_in_block;
         Size_send_U = Size_send_U + rows_in_block*cols_in_block;
         cols_in_buffer_U_my_initial = cols_in_buffer_U_my_initial + cols_in_block;
      }
      curr_col_loc = curr_col_loc + nblk;
      rows_in_block = rows_in_block + ratio*nblk;
   }
   rows_in_buffer_U_my_initial = rows_in_block - ratio*nblk;
   *Buf_pos = (double complex)cols_in_buffer_U_my_initial;
   Buf_pos = Buf_pos + 1;
   *Buf_pos = (double complex)rows_in_buffer_U_my_initial;
   Size_send_U = Size_send_U + 2;

   proc_col_min = np_cols;
   for(i = 0; i < ratio; i++)
   {
      from_where_to_receive_U = (my_pcol + my_prow + i*np_rows)%np_cols;
      if(from_where_to_receive_U < proc_col_min)
         proc_col_min = from_where_to_receive_U;
   }

   Size_receive_U = 0;
   cols_in_buffer_U = 0;
   rows_in_buffer_U = 0;
   for(i = 0; i < ratio; i++)
   {
      where_to_send_U = (my_pcol - my_prow - i*np_rows + np_cols)%np_cols;
      from_where_to_receive_U = (my_pcol + my_prow + i*np_rows)%np_cols;

      if(ratio != 1)
      {
         if(where_to_send_U != my_pcol)
         {
            MPI_Sendrecv(Buf_to_send_U, (int) Size_send_U, MPI_DOUBLE_COMPLEX, (int) where_to_send_U, 0, Buf_U, (int) Size_U_stored, MPI_DOUBLE_COMPLEX, (int) from_where_to_receive_U, 0, row_comm, &status);
            MPI_Get_count(&status, MPI_DOUBLE_COMPLEX, &Size_receive_U_nowMPI);
            Size_receive_U_now = (int) Size_receive_U_nowMPI;
            Size_receive_U = Size_receive_U + Size_receive_U_now - 2;

            cols_in_buffer_U_now = Buf_U[Size_receive_U_now - 2];
            cols_in_buffer_U = cols_in_buffer_U + cols_in_buffer_U_now;
            rows_in_buffer_U_now = Buf_U[Size_receive_U_now - 1];

            if(rows_in_buffer_U < rows_in_buffer_U_now)
               rows_in_buffer_U = rows_in_buffer_U_now;

            intNumber = from_where_to_receive_U/np_rows;
            if(proc_col_min >= my_prow)
               CopyTo = &Buf_to_receive_U[nblk*nblk*intNumber*(intNumber + 1)/2];
            else
               if(from_where_to_receive_U < my_prow)
                  CopyTo = &Buf_to_receive_U[nblk*nblk*ratio*(ratio - 1)/2];
               else
                  CopyTo = &Buf_to_receive_U[nblk*nblk*intNumber*(intNumber - 1)/2];
            CopyFrom = Buf_U;
         }
         else
         {
            cols_in_buffer_U_now = cols_in_buffer_U_my_initial;
            cols_in_buffer_U = cols_in_buffer_U + cols_in_buffer_U_now;

            rows_in_buffer_U_now = rows_in_buffer_U_my_initial;
            if(rows_in_buffer_U < rows_in_buffer_U_now)
               rows_in_buffer_U = rows_in_buffer_U_now;

            intNumber = my_pcol/np_rows;
            if(proc_col_min >= my_prow)
               CopyTo = &Buf_to_receive_U[nblk*nblk*intNumber*(intNumber + 1)/2];
            else
               if(my_pcol < my_prow)
                  CopyTo = &Buf_to_receive_U[nblk*nblk*ratio*(ratio - 1)/2];
               else
                  CopyTo = &Buf_to_receive_U[nblk*nblk*intNumber*(intNumber - 1)/2];
            CopyFrom = Buf_to_send_U;
            Size_receive_U = Size_receive_U + Size_send_U - 2;
         }

         intNumber = ceil((double complex)cols_in_buffer_U_now/(double complex)nblk);
         if(from_where_to_receive_U >= my_prow)
            rows_in_block = ceil(((double complex)(from_where_to_receive_U + 1) - (double complex)my_prow)/(double complex)np_rows)*nblk;
         else
            rows_in_block = ratio*nblk;
         for(j = 0; j < intNumber; j++)
         {
            if((j+1)*nblk < cols_in_buffer_U_now)
               cols_in_block = nblk;
            else
               cols_in_block = cols_in_buffer_U_now - j*nblk;

            zlacpy("A", &rows_in_block, &cols_in_block, CopyFrom, &rows_in_block, CopyTo, &rows_in_block);

            CopyFrom = CopyFrom + rows_in_block*cols_in_block;
            CopyTo = CopyTo + ratio*rows_in_block*nblk + nblk*nblk*ratio*(ratio-1)/2;
            rows_in_block = rows_in_block + ratio*nblk;
            if(rows_in_block > rows_in_buffer_U_now)
               rows_in_block = rows_in_buffer_U_now;
         }
      }
      else
      {
         if(my_prow > 0)
         {
            MPI_Sendrecv(Buf_to_send_U, (int) Size_send_U, MPI_DOUBLE_COMPLEX, (int) where_to_send_U, 0, Buf_to_receive_U, (int) Size_U_stored, MPI_DOUBLE_COMPLEX, (int) from_where_to_receive_U, 0, row_comm, &status);
            MPI_Get_count(&status, MPI_DOUBLE_COMPLEX, &Size_receive_UMPI);
            Size_receive_U = (int) Size_receive_UMPI;

            cols_in_buffer_U = (int)Buf_to_receive_U[Size_receive_U-2];
            rows_in_buffer_U = (int)Buf_to_receive_U[Size_receive_U-1];
         }
         else
         {
            Size_receive_U = Size_send_U;
            rows_in_buffer_U = rows_in_buffer_U_my_initial;
            cols_in_buffer_U = cols_in_buffer_U_my_initial;
         }
      }
   }
   if(ratio != 1)
   {
      Buf_to_receive_U[Size_receive_U] = cols_in_buffer_U;
      Buf_to_receive_U[Size_receive_U + 1] = rows_in_buffer_U;
      Size_receive_U = Size_receive_U + 2;
   }

   if(my_pcol > 0)
   {
      where_to_send_B = (my_prow - my_pcol + np_cols)%np_rows;
      from_where_to_receive_B = (my_pcol + my_prow)%np_rows;

      if(where_to_send_B != my_prow)
      {

         zlacpy("A", &na_rows, &nb_cols, B, &na_rows, Buf_to_send_B, &na_rows);
         MPI_Sendrecv(Buf_to_send_B, (int) (nb_cols*na_rows), MPI_DOUBLE_COMPLEX, (int) where_to_send_B, 0, Buf_to_receive_B, (int) (nb_cols*Buf_rows), MPI_DOUBLE_COMPLEX, (int) from_where_to_receive_B, 0, col_comm, &status);
         MPI_Get_count(&status, MPI_DOUBLE_COMPLEX, &Size_receive_BMPI);
         Size_receive_B = (int) Size_receive_BMPI;
         Size_receive_B = Size_receive_B/nb_cols;

      }
      else
      {
         zlacpy("A", &na_rows, &nb_cols, B, &na_rows, Buf_to_receive_B, &na_rows);
         Size_receive_B = na_rows;
      }
   }
   else
   {
      zlacpy("A", &na_rows, &nb_cols, B, &na_rows, Buf_to_receive_B, &na_rows);
      Size_receive_B = na_rows;
   }

   where_to_send_U = (my_pcol - 1 + np_cols)%np_cols;
   from_where_to_receive_U = (my_pcol + 1)%np_cols;
   where_to_send_B = (my_prow - 1 + np_rows)%np_rows;
   from_where_to_receive_B = (my_prow + 1)%np_rows;

   for(i = 1; i < np_rows; i++)
   {

      double_ptr = Buf_to_send_U;
      Buf_to_send_U = Buf_to_receive_U;
      Buf_to_receive_U = double_ptr;

      double_ptr = Buf_to_send_B;
      Buf_to_send_B = Buf_to_receive_B;
      Buf_to_receive_B = double_ptr;

      Size_send_U = Size_receive_U;
      Size_send_B = Size_receive_B;

      MPI_Isend(Buf_to_send_U, (int) Size_send_U, MPI_DOUBLE_COMPLEX, (int) where_to_send_U, 0, row_comm, &request_U_Send);
      MPI_Irecv(Buf_to_receive_U, (int) (ratio*Size_U_stored), MPI_DOUBLE_COMPLEX, (int) from_where_to_receive_U, 0, row_comm, &request_U_Recv);

      MPI_Isend(Buf_to_send_B, (int) (Size_send_B*nb_cols), MPI_DOUBLE_COMPLEX, (int) where_to_send_B, 0, col_comm, &request_B_Send);
      MPI_Irecv(Buf_to_receive_B, (int) (Buf_rows*nb_cols), MPI_DOUBLE_COMPLEX, (int) from_where_to_receive_B, 0, col_comm, &request_B_Recv);

      cols_in_buffer_U = (int)Buf_to_send_U[Size_receive_U-2];
      rows_in_buffer_U = (int)Buf_to_send_U[Size_receive_U-1];

      proc_col_min = np_cols;
      for(j = 0; j < ratio; j++)
      {
         col_of_origin_U = (my_pcol + my_prow + i - 1 + j*np_rows)%np_cols;
         if(col_of_origin_U < proc_col_min)
            proc_col_min = col_of_origin_U;
      }
      col_of_origin_U = proc_col_min;

      num_of_blocks_in_U_buffer = ceil((double complex)cols_in_buffer_U/(double complex)nblk);

      if (col_of_origin_U >= my_prow)
         B_local_start = Buf_to_send_B;
      else
         B_local_start = Buf_to_send_B + nblk;

      U_local_start = Buf_to_send_U;

      for(j = 0; j < num_of_blocks_in_U_buffer; j++)
      {
         curr_rows = (j+1)*nblk;
         if (curr_rows > rows_in_buffer_U)
            curr_rows = rows_in_buffer_U;

         if((j+1)*nblk <= cols_in_buffer_U)
            b_rows_mult = nblk;
         else
            b_rows_mult = cols_in_buffer_U - j*nblk;

         zgemm_("N", "N", &curr_rows, &nb_cols, &b_rows_mult, &done, U_local_start, &curr_rows, B_local_start, &Size_receive_B, &done, Res, &na_rows);

         U_local_start = U_local_start + nblk*curr_rows;
         B_local_start = B_local_start + nblk;
      }

      MPI_Wait(&request_U_Send, &status);
      MPI_Wait(&request_U_Recv, &status);
      MPI_Get_count(&status, MPI_DOUBLE_COMPLEX, &Size_receive_UMPI);
      Size_receive_U = (int) Size_receive_UMPI;

      MPI_Wait(&request_B_Send, &status);
      MPI_Wait(&request_B_Recv, &status);
      MPI_Get_count(&status, MPI_DOUBLE_COMPLEX, &Size_receive_BMPI);
      Size_receive_B = (int) Size_receive_BMPI;
      Size_receive_B = (int) Size_receive_B / nb_cols;

   }

   cols_in_buffer_U = (int)Buf_to_receive_U[Size_receive_U-2];
   rows_in_buffer_U = (int)Buf_to_receive_U[Size_receive_U-1];

   proc_col_min = np_cols;
   for(j = 0; j < ratio; j++)
   {
      col_of_origin_U = (my_pcol + my_prow + np_rows - 1 + j*np_rows)%np_cols;
      if(col_of_origin_U < proc_col_min)
         proc_col_min = col_of_origin_U;
   }
   col_of_origin_U = proc_col_min;

   num_of_blocks_in_U_buffer = ceil((double complex)cols_in_buffer_U/(double complex)nblk);

   if (col_of_origin_U >= my_prow)
      B_local_start = Buf_to_receive_B;
   else
      B_local_start = Buf_to_receive_B + nblk;

   U_local_start = Buf_to_receive_U;

   for(j = 0; j < num_of_blocks_in_U_buffer; j++)
   {
      curr_rows = (j+1)*nblk;
      if (curr_rows > rows_in_buffer_U)
         curr_rows = rows_in_buffer_U;

      if((j+1)*nblk <= cols_in_buffer_U)
         b_rows_mult = nblk;
      else
         b_rows_mult = cols_in_buffer_U - j*nblk;

      zgemm_("N", "N", &curr_rows, &nb_cols, &b_rows_mult, &done, U_local_start, &curr_rows, B_local_start, &Size_receive_B, &done, Res, &na_rows);

      U_local_start = U_local_start + nblk*curr_rows;
      B_local_start = B_local_start + nblk;
   }

   free(Buf_to_send_U);
   free(Buf_to_receive_U);
   free(Buf_to_send_B);
   free(Buf_to_receive_B);
   if(ratio != 1)
      free(Buf_U);
}

void cannons_triang_rectangular_c_dc(double complex* U, double complex* B, int local_rowsCast, int local_colsCast,
                                    int* u_desc, int* b_desc, double complex *Res, int row_comm, int col_comm)
{
  int local_rows, local_cols;

  local_rows = (int) local_rowsCast;
  local_cols = (int) local_colsCast;

  MPI_Comm c_row_comm = MPI_Comm_f2c(row_comm);
  MPI_Comm c_col_comm = MPI_Comm_f2c(col_comm);

  int my_prow, my_pcol, np_rows, np_cols;
  int my_prowMPI, my_pcolMPI, np_rowsMPI, np_colsMPI;

  MPI_Comm_rank(c_row_comm, &my_prowMPI);
  MPI_Comm_size(c_row_comm, &np_rowsMPI);
  MPI_Comm_rank(c_col_comm, &my_pcolMPI);
  MPI_Comm_size(c_col_comm, &np_colsMPI);

  my_prow = (int) my_prowMPI;
  my_pcol = (int) my_pcolMPI;
  np_rows = (int) np_rowsMPI;
  np_cols = (int) np_colsMPI;

  cannons_triang_rectangular_dc(U, B, np_rows, np_cols, my_prow, my_pcol, u_desc, b_desc, Res, c_col_comm, c_row_comm);
}
void cannons_reduction_c_dc(double complex* A, double complex* U, int local_rowsCast, int local_colsCasr, int* a_desc,
                            double complex *Res, int ToStore, int row_comm, int col_comm);
void cannons_triang_rectangular_c_dc(double complex* U, double complex* B, int local_rowsCast, int local_colsCast,
                                    int* u_desc, int* b_desc, double complex *Res, int row_comm, int col_comm);

void dlacpy(char*, int*, int*, double*, int*, double*, int*);
void dgemm_(char*, char*, int*, int*, int*, double*, double*, int*, double*, int*, double*, double*, int*);

void slacpy(char*, int*, int*, float*, int*, float*, int*);
void sgemm_(char*, char*, int*, int*, int*, float*, float*, int*, float*, int*, float*, float*, int*);

void zlacpy(char*, int*, int*, double complex*, int*, double complex*, int*);
void zgemm_(char*, char*, int*, int*, int*, double complex*, double complex*, int*, double complex*, int*, double complex*, double complex*, int*);

void clacpy(char*, int*, int*, float complex*, int*, float complex*, int*);
void cgemm_(char*, char*, int*, int*, int*, float complex*, float complex*, int*, float complex*, int*, float complex*, float complex*, int*);
int numroc(int*, int*, int*, int*, int*);

void pdlacpy(char*, int*, int*, double*, int*, int*, int*, double*, int*, int*, int*);
void pdtran(int*, int*, double*, double*, int*, int*, int*, double*, double*, int*, int*, int*);

void pslacpy(char*, int*, int*, float*, int*, int*, int*, float*, int*, int*, int*);
void pstran(int*, int*, float*, float*, int*, int*, int*, float*, float*, int*, int*, int*);

void pzlacpy(char*, int*, int*, double complex*, int*, int*, int*, double complex*, int*, int*, int*);
void pztranc(int*, int*, double complex*, double complex*, int*, int*, int*, double complex*, double complex*, int*, int*, int*);

void pclacpy(char*, int*, int*, float complex*, int*, int*, int*, float complex*, int*, int*, int*);
void pctranc(int*, int*, float complex*, float complex*, int*, int*, int*, float complex*, float complex*, int*, int*, int*);

void cannons_reduction_fc(float complex* A, float complex* U, int np_rows, int np_cols, int my_prow, int my_pcol,
                         int* a_desc, float complex *Res, int ToStore, MPI_Comm row_comm, MPI_Comm col_comm)
{
   int na, nblk, i, j, Size_send_A, Size_receive_A, Size_send_U, Size_receive_U, Buf_rows, Buf_cols, where_to_send_A, from_where_to_receive_A, where_to_send_U, from_where_to_receive_U, last_proc_row, last_proc_col, cols_in_buffer_A, rows_in_buffer_A, intNumber;
   float complex *Buf_to_send_A, *Buf_to_receive_A, *Buf_to_send_U, *Buf_to_receive_U, *data_ptr, *Buf_A, *Buf_pos, *U_local_start, *Res_ptr, *M, *M_T, *A_local_start, *U_local_start_curr, *U_stored, *CopyTo, *CopyFrom, *U_to_calc;
   int ratio, num_of_iters, cols_in_buffer, rows_in_block, rows_in_buffer, curr_col_loc, cols_in_block, curr_col_glob, curr_row_loc, Size_receive_A_now, Nb, owner, cols_in_buffer_A_now;
   int Size_receive_A_nowMPI, Size_receive_AMPI, Size_receive_UMPI;

   int row_of_origin_U, rows_in_block_U, num_of_blocks_in_U_buffer, k, startPos, cols_in_buffer_U, rows_in_buffer_U, col_of_origin_A, curr_row_loc_res, curr_row_loc_A, curr_col_glob_res;
   int curr_col_loc_res, curr_col_loc_buf, proc_row_curr, curr_col_loc_U, A_local_index, LDA_A, LDA_A_new, index_row_A_for_LDA, ii, rows_in_block_U_curr, width, row_origin_U, rows_in_block_A, cols_in_buffer_A_my_initial, rows_in_buffer_A_my_initial, proc_col_min;
   int *SizesU;
   int Size_U_skewed, Size_U_stored, Curr_pos_in_U_stored, rows_in_buffer_A_now;
   float complex done = 1.0;
   float complex dzero = 0.0;
   int one = 1;
   int zero = 0;
   int na_rows, na_cols;

   MPI_Status status;
   MPI_Request request_A_Recv;
   MPI_Request request_A_Send;
   MPI_Request request_U_Recv;
   MPI_Request request_U_Send;

   na = a_desc[2];
   nblk = a_desc[4];
   na_rows = numroc(&na, &nblk, &my_prow, &zero, &np_rows);
   na_cols = numroc(&na, &nblk, &my_pcol, &zero, &np_cols);

   if(ToStore > (np_rows -1))
      if((my_prow == 0)&&(my_pcol == 0))
         printf("Buffering level is larger than (np_rows-1) !!!\n");
   if((my_prow == 0)&&(my_pcol == 0))
         printf("Buffering level = %d\n", ToStore);

   if (np_cols%np_rows != 0)
   {
      if((my_prow == 0)&& (my_pcol ==0))
         printf("!!!!! np_cols must be a multiple of np_rows!!!!! I do nothing! \n");
      return;
   }
   if (np_cols < np_rows != 0)
   {
      if((my_prow == 0)&& (my_pcol ==0))
         printf("np_cols < np_rows \n");
      return;
   }

   ratio = np_cols/np_rows;
   last_proc_row = ((na-1)/nblk) % np_rows;
   last_proc_col = ((na-1)/nblk) % np_cols;

   if(na%nblk == 0)
      if(my_pcol <= last_proc_col)
         Buf_cols = na_cols;
      else
         Buf_cols = na_cols + nblk;
   else
      if(my_pcol < last_proc_col)
         Buf_cols = na_cols;
      else if(my_pcol > last_proc_col)
         Buf_cols = na_cols + nblk;
      else
         Buf_cols = na_cols + nblk - na_cols%nblk;

  if(na%nblk == 0)
      if(my_prow <= last_proc_row)
         Buf_rows = na_rows;
      else
         Buf_rows = na_rows + nblk;
   else
      if(my_prow < last_proc_row)
         Buf_rows = na_rows;
      else if(my_prow > last_proc_row)
         Buf_rows = na_rows + nblk;
      else
         Buf_rows = na_rows + nblk - na_rows%nblk;

   intNumber = ceil((float complex)na/(float complex)(np_cols*nblk));
   Size_U_stored = ratio*nblk*nblk*intNumber*(intNumber+1)/2 + 2;

   U_stored = malloc((Size_U_stored*(ToStore+1))*sizeof(float complex));
   SizesU = malloc(ToStore*sizeof(int));
   Buf_to_send_A = malloc(ratio*Buf_cols*Buf_rows*sizeof(float complex));
   Buf_to_receive_A = malloc(ratio*Buf_cols*Buf_rows*sizeof(float complex));
   Buf_to_send_U = malloc(Size_U_stored*sizeof(float complex));
   Buf_to_receive_U = malloc(Size_U_stored*sizeof(float complex));
   if(ratio != 1)
      Buf_A = malloc(Buf_cols*Buf_rows*sizeof(float complex));
   M = malloc(na_rows*na_cols*sizeof(float complex));
   M_T = malloc(na_rows*na_cols*sizeof(float complex));
   for(i = 0; i < na_rows*na_cols; i++)
      M[i] = 0;

   if(ratio != 1)
      clacpy("A", &na_rows, &na_cols, A, &na_rows, Buf_to_send_A, &na_rows);
   Size_receive_A = 0;

   for(i = 0; i < ratio; i++)
   {
      where_to_send_A = (my_pcol - my_prow - i*np_rows + np_cols)%np_cols;
      from_where_to_receive_A = (my_pcol + my_prow + i*np_rows)%np_cols;

      if(ratio != 1)
      {
         if(where_to_send_A != my_pcol)
         {
           MPI_Sendrecv(Buf_to_send_A, (int) (na_cols*na_rows), MPI_COMPLEX,(int) where_to_send_A, (int) zero, Buf_A, (int) (na_rows*Buf_cols), MPI_COMPLEX, (int) from_where_to_receive_A, (int) zero, row_comm, &status);
           MPI_Get_count(&status, MPI_COMPLEX, &Size_receive_A_nowMPI);
           Size_receive_A_now = (int) Size_receive_A_nowMPI;
           Size_receive_A_now = Size_receive_A_now/na_rows;
         }
         else
            Size_receive_A_now = na_cols;
         Size_receive_A = Size_receive_A + Size_receive_A_now;

         intNumber = from_where_to_receive_A/np_rows;

         CopyTo = &Buf_to_receive_A[intNumber*na_rows*nblk];
         if(where_to_send_A != my_pcol)
            CopyFrom = Buf_A;
         else
            CopyFrom = A;

         intNumber = ceil((float complex)Size_receive_A_now/(float complex)nblk);
         for(j = 0; j < intNumber; j++)
         {
            width = nblk;
            if(nblk*(j+1) > Size_receive_A_now)
               width = Size_receive_A_now - nblk*j;
            clacpy("A", &na_rows, &width, CopyFrom, &na_rows, CopyTo, &na_rows);
            CopyTo = CopyTo + na_rows*nblk*ratio;
            CopyFrom = CopyFrom + na_rows*nblk;
         }
      }
      else
         if(my_prow > 0)
         {
            clacpy("A", &na_rows, &na_cols, A, &na_rows, Buf_to_send_A, &na_rows);
            MPI_Sendrecv(Buf_to_send_A, (int) (na_cols*na_rows), MPI_COMPLEX, (int) where_to_send_A, (int) zero, Buf_to_receive_A, (int) (na_rows*Buf_cols), MPI_COMPLEX, (int) from_where_to_receive_A, (int) zero, row_comm, &status);
            MPI_Get_count(&status, MPI_COMPLEX, &Size_receive_AMPI);
            Size_receive_A = (int) Size_receive_AMPI;
     Size_receive_A = Size_receive_A/na_rows;
         }
         else
         {
            clacpy("A", &na_rows, &na_cols, A, &na_rows, Buf_to_receive_A, &na_rows);
            Size_receive_A = na_cols;
         }
   }

   num_of_iters = ceil((float complex)na_cols/(float complex)nblk);

   where_to_send_U = (my_prow - my_pcol + np_cols)%np_rows;
   from_where_to_receive_U = (my_pcol + my_prow)%np_rows;

   if(where_to_send_U == my_prow)
      Buf_pos = Buf_to_receive_U;
   else
      Buf_pos = Buf_to_send_U;

   if(my_pcol >= my_prow)
      curr_col_loc = 0;
   else
      curr_col_loc = 1;

   num_of_iters = num_of_iters - curr_col_loc;
   curr_col_loc = curr_col_loc*nblk;

   if(my_pcol >= my_prow )
      rows_in_block = ceil(((float complex)(my_pcol + 1) - (float complex)my_prow)/(float complex)np_rows)*nblk;
   else
      rows_in_block = ratio*nblk;

   Size_send_U = 0;
   for(i = 0; i < num_of_iters; i++)
   {
      if(rows_in_block > na_rows)
         rows_in_block = na_rows;

      if ((na_cols - curr_col_loc) < nblk)
         cols_in_block = na_cols - curr_col_loc;
      else
         cols_in_block = nblk;

      if((rows_in_block > 0)&&(cols_in_block > 0))
      {
         data_ptr = &U[curr_col_loc*na_rows];
         clacpy("A", &rows_in_block, &cols_in_block, data_ptr, &na_rows, Buf_pos, &rows_in_block);
         Buf_pos = Buf_pos + rows_in_block*cols_in_block;
         Size_send_U = Size_send_U + rows_in_block*cols_in_block;
      }
      curr_col_loc = curr_col_loc + nblk;
      rows_in_block = rows_in_block + ratio*nblk;
   }
   rows_in_buffer = rows_in_block - ratio*nblk;
   *Buf_pos = (float complex)rows_in_buffer;
   Size_send_U = Size_send_U + 1;

   if(where_to_send_U != my_prow)
   {

      MPI_Sendrecv(Buf_to_send_U, (int) Size_send_U, MPI_COMPLEX, (int) where_to_send_U, (int) zero, Buf_to_receive_U, (int) (Buf_rows*na_cols), MPI_COMPLEX, (int) from_where_to_receive_U, (int) zero, col_comm, &status);
      MPI_Get_count(&status, MPI_COMPLEX, &Size_receive_UMPI);
      Size_receive_U = (int) Size_receive_UMPI;
   }
   else
      Size_receive_U = Size_send_U;

   for(i = 0; i < Size_receive_U; i++)
      U_stored[i] = Buf_to_receive_U[i];
   Size_U_skewed = Size_receive_U;
   Curr_pos_in_U_stored = Size_U_skewed;

   where_to_send_A = (my_pcol - 1 + np_cols)%np_cols;
   from_where_to_receive_A = (my_pcol + 1)%np_cols;
   where_to_send_U = (my_prow - 1 + np_rows)%np_rows;
   from_where_to_receive_U = (my_prow + 1)%np_rows;

   for(j = 1; j < np_rows; j++)
   {

      data_ptr = Buf_to_send_A;
      Buf_to_send_A = Buf_to_receive_A;
      Buf_to_receive_A = data_ptr;

      data_ptr = Buf_to_send_U;
      Buf_to_send_U = Buf_to_receive_U;
      Buf_to_receive_U = data_ptr;

      Size_send_A = Size_receive_A;
      MPI_Isend(Buf_to_send_A, (int) (Size_send_A*na_rows), MPI_COMPLEX, (int) where_to_send_A, (int) zero, row_comm, &request_A_Send);
      MPI_Irecv(Buf_to_receive_A, (int) (Buf_cols*na_rows*ratio), MPI_COMPLEX, (int) from_where_to_receive_A, (int) zero, row_comm, &request_A_Recv);

      Size_send_U = Size_receive_U;
      MPI_Isend(Buf_to_send_U, (int) Size_send_U, MPI_COMPLEX, (int) where_to_send_U, (int) zero, col_comm, &request_U_Send);
      MPI_Irecv(Buf_to_receive_U, (int) (Buf_rows*na_cols), MPI_COMPLEX, (int) from_where_to_receive_U, (int) zero, col_comm, &request_U_Recv);

      rows_in_buffer = (int)Buf_to_send_U[Size_receive_U-1];
      row_origin_U = (my_pcol + my_prow + np_cols + j - 1)%np_rows;

      if((my_pcol >= my_prow)&&(my_pcol >= row_origin_U))
      {
         cols_in_buffer = na_cols;
         curr_col_loc_res = 0;
         curr_col_loc_buf = 0;
      }
      if((my_pcol < my_prow)&&(my_pcol < row_origin_U))
      {
         cols_in_buffer = na_cols - nblk;
         curr_col_loc_res = nblk;
         curr_col_loc_buf = 0;
      }
      if((my_pcol >= my_prow)&&(my_pcol < row_origin_U))
      {
         cols_in_buffer = na_cols - nblk;
         curr_col_loc_res = nblk;
         curr_col_loc_buf = 0;
      }
      if((my_pcol < my_prow)&&(my_pcol >= row_origin_U))
      {
         cols_in_buffer = na_cols;
         curr_col_loc_res = nblk;
         curr_col_loc_buf = nblk;
      }

      num_of_blocks_in_U_buffer = ceil(((float complex)cols_in_buffer - (float complex)curr_col_loc_buf)/(float complex)nblk);

      startPos = (curr_col_loc_buf + nblk)*curr_col_loc_buf/2;
      U_local_start = &Buf_to_send_U[startPos];
      Res_ptr = &M[curr_col_loc_res*na_rows];

      for (i = 0; i < num_of_blocks_in_U_buffer; i++)
      {
         curr_col_glob = (curr_col_loc_res/nblk)*nblk*np_cols + my_pcol*nblk;
         proc_row_curr = (curr_col_glob/nblk)%np_rows;
         rows_in_block_A = (curr_col_glob/(nblk*np_rows))*nblk;
         if(my_prow <= proc_row_curr)
            rows_in_block_A = rows_in_block_A + nblk;

         if(rows_in_block_A > na_rows)
            rows_in_block_A = na_rows;

         if((curr_col_loc_buf + nblk) <= cols_in_buffer)
            cols_in_block = nblk;
         else
            cols_in_block = cols_in_buffer - curr_col_loc_buf;

         rows_in_block_U = (curr_col_glob/(nblk*np_rows))*nblk;
         if(proc_row_curr >= row_origin_U)
            rows_in_block_U = rows_in_block_U + nblk;

         if(rows_in_block_U > rows_in_buffer)
            rows_in_block_U = rows_in_buffer;

         if ((rows_in_block_A > 0)&&(cols_in_block > 0))
            if (j == 1) {
               cgemm_("N", "N", &rows_in_block_A, &cols_in_block, &rows_in_block_U, &done, Buf_to_send_A, &na_rows, U_local_start, &rows_in_block_U, &dzero, Res_ptr, &na_rows);
     }
            else {
               cgemm_("N", "N", &rows_in_block_A, &cols_in_block, &rows_in_block_U, &done, Buf_to_send_A, &na_rows, U_local_start, &rows_in_block_U, &done, Res_ptr, &na_rows);
     }

         U_local_start = U_local_start + rows_in_block_U*cols_in_block;
         curr_col_loc_res = curr_col_loc_res + nblk;
         Res_ptr = &M[curr_col_loc_res*na_rows];
         curr_col_loc_buf = curr_col_loc_buf + nblk;
      }

      MPI_Wait(&request_A_Send, &status);
      MPI_Wait(&request_A_Recv, &status);

      MPI_Get_count(&status, MPI_COMPLEX, &Size_receive_AMPI);
      Size_receive_A = (int) Size_receive_AMPI;
      Size_receive_A = Size_receive_A / na_rows;

      MPI_Wait(&request_U_Send, &status);
      MPI_Wait(&request_U_Recv, &status);
      MPI_Get_count(&status, MPI_COMPLEX, &Size_receive_UMPI);
      Size_receive_U = (int) Size_receive_UMPI;

      if(j <= ToStore)
      {
         for(k = 0; k < Size_receive_U; k++)
            U_stored[Curr_pos_in_U_stored + k] = Buf_to_receive_U[k];
         Curr_pos_in_U_stored = Curr_pos_in_U_stored + Size_receive_U;
         SizesU[j-1] = Size_receive_U;
      }
   }

   rows_in_buffer = (int)Buf_to_receive_U[Size_receive_U-1];
   row_origin_U = (my_pcol + my_prow + np_cols + np_rows -1)%np_rows;

   if((my_pcol >= my_prow)&&(my_pcol >= row_origin_U))
   {
      cols_in_buffer = na_cols;
      curr_col_loc_res = 0;
      curr_col_loc_buf = 0;
   }
   if((my_pcol < my_prow)&&(my_pcol < row_origin_U))
   {
      cols_in_buffer = na_cols - nblk;
      curr_col_loc_res = nblk;
      curr_col_loc_buf = 0;
   }
   if((my_pcol >= my_prow)&&(my_pcol < row_origin_U))
   {
      cols_in_buffer = na_cols - nblk;
      curr_col_loc_res = nblk;
      curr_col_loc_buf = 0;
   }
   if((my_pcol < my_prow)&&(my_pcol >= row_origin_U))
   {
      cols_in_buffer = na_cols;
      curr_col_loc_res = nblk;
      curr_col_loc_buf = nblk;
   }

   num_of_blocks_in_U_buffer = ceil(((float complex)cols_in_buffer - (float complex)curr_col_loc_buf)/(float complex)nblk);

   startPos = (curr_col_loc_buf + nblk)*curr_col_loc_buf/2;
   U_local_start = &Buf_to_receive_U[startPos];
   Res_ptr = &M[curr_col_loc_res*na_rows];

   for (i = 0; i < num_of_blocks_in_U_buffer; i++)
   {
      curr_col_glob = (curr_col_loc_res/nblk)*nblk*np_cols + my_pcol*nblk;
      proc_row_curr = (curr_col_glob/nblk)%np_rows;
      rows_in_block_A = (curr_col_glob/(nblk*np_rows))*nblk;
      if(my_prow <= proc_row_curr)
         rows_in_block_A = rows_in_block_A + nblk;

      if(rows_in_block_A > na_rows)
         rows_in_block_A = na_rows;

      if((curr_col_loc_buf + nblk) <= cols_in_buffer)
         cols_in_block = nblk;
      else
         cols_in_block = cols_in_buffer - curr_col_loc_buf;

      rows_in_block_U = (curr_col_glob/(nblk*np_rows))*nblk;
      if(proc_row_curr >= row_origin_U)
         rows_in_block_U = rows_in_block_U + nblk;

      if(rows_in_block_U > rows_in_buffer)
         rows_in_block_U = rows_in_buffer;

      if ((rows_in_block_A > 0)&&(cols_in_block > 0))
         if (j == 1) {
            cgemm_("N", "N", &rows_in_block_A, &cols_in_block, &rows_in_block_U, &done, Buf_to_receive_A, &na_rows, U_local_start, &rows_in_block_U, &dzero, Res_ptr, &na_rows);
  }
         else {
            cgemm_("N", "N", &rows_in_block_A, &cols_in_block, &rows_in_block_U, &done, Buf_to_receive_A, &na_rows, U_local_start, &rows_in_block_U, &done, Res_ptr, &na_rows);
         }

      U_local_start = U_local_start + rows_in_block_U*cols_in_block;
      curr_col_loc_res = curr_col_loc_res + nblk;
      Res_ptr = &M[curr_col_loc_res*na_rows];
      curr_col_loc_buf = curr_col_loc_buf + nblk;
   }

   pctranc(&na, &na, &done, M, &one, &one, a_desc, &dzero, M_T, &one, &one, a_desc);

   if((ratio != 1)||(my_prow != 0))
      Buf_pos = Buf_to_send_A;
   else
      Buf_pos = Buf_to_receive_A;

   num_of_iters = ceil((float complex)na_cols/(float complex)nblk);

   cols_in_buffer_A_my_initial = 0;
   Size_send_A = 0;

   if(my_pcol <= my_prow)
   {
      curr_row_loc = 0;
      rows_in_buffer_A_my_initial = na_rows;
   }
   else
   {
      curr_row_loc = ceil((float complex)(((float complex)my_pcol - (float complex)my_prow)/(float complex)np_rows))*nblk;
      rows_in_buffer_A_my_initial = na_rows - curr_row_loc;
   }

   for(i = 0; i < num_of_iters; i++)
   {
      curr_col_loc = i*nblk;
      rows_in_block = na_rows - curr_row_loc;

      if ((na_cols - curr_col_loc) < nblk)
         cols_in_block = na_cols - curr_col_loc;
      else
         cols_in_block = nblk;

      if((rows_in_block > 0)&&(cols_in_block > 0))
      {
         A_local_start = &M_T[curr_col_loc*na_rows + curr_row_loc];
         clacpy("A", &rows_in_block, &cols_in_block, A_local_start, &na_rows, Buf_pos, &rows_in_block);
         Buf_pos = Buf_pos + rows_in_block*cols_in_block;
         Size_send_A = Size_send_A + rows_in_block*cols_in_block;
         cols_in_buffer_A_my_initial = cols_in_buffer_A_my_initial + cols_in_block;
      }
      curr_row_loc = curr_row_loc + ratio*nblk;
   }
   *Buf_pos = (float complex)cols_in_buffer_A_my_initial;
   Size_send_A = Size_send_A + 1;

   proc_col_min = np_cols;
   for(i = 0; i < ratio; i++)
   {
      from_where_to_receive_A = (my_pcol + my_prow + i*np_rows)%np_cols;
      if(from_where_to_receive_A < proc_col_min)
         proc_col_min = from_where_to_receive_A;
   }

   Size_receive_A = 0;
   cols_in_buffer_A = 0;
   rows_in_buffer_A = 0;
   for(i = 0; i < ratio; i++)
   {
      where_to_send_A = (my_pcol - my_prow - i*np_rows + np_cols)%np_cols;
      from_where_to_receive_A = (my_pcol + my_prow + i*np_rows)%np_cols;

      if(ratio != 1)
      {
         if(where_to_send_A != my_pcol)
         {
            MPI_Sendrecv(Buf_to_send_A, (int) Size_send_A, MPI_COMPLEX, (int) where_to_send_A, (int) zero, Buf_A, (int) Size_U_stored, MPI_COMPLEX, (int) from_where_to_receive_A, (int) zero, row_comm, &status);
            MPI_Get_count(&status, MPI_COMPLEX, &Size_receive_A_nowMPI);
            Size_receive_A_now = (int) Size_receive_A_nowMPI;

            Size_receive_A = Size_receive_A + Size_receive_A_now - 1;

            cols_in_buffer_A_now = Buf_A[Size_receive_A_now-1];
            cols_in_buffer_A = cols_in_buffer_A + cols_in_buffer_A_now;

            if(from_where_to_receive_A <= my_prow)
            {
               rows_in_buffer_A_now = na_rows;
            }
            else
            {
               rows_in_buffer_A_now = na_rows - ceil((float complex)(((float complex)from_where_to_receive_A - (float complex)my_prow)/(float complex)np_rows))*nblk;
            }
            if(rows_in_buffer_A < rows_in_buffer_A_now)
               rows_in_buffer_A = rows_in_buffer_A_now;

            intNumber = from_where_to_receive_A/np_rows;
            if(proc_col_min <= my_prow)
               CopyTo = &Buf_to_receive_A[nblk*(na_rows*intNumber - nblk*(intNumber-1)*intNumber/2)];
            else
               CopyTo = &Buf_to_receive_A[nblk*(na_rows*intNumber - nblk*intNumber*(intNumber+1)/2)];
            CopyFrom = Buf_A;
         }
         else
         {
            cols_in_buffer_A_now = cols_in_buffer_A_my_initial;
            cols_in_buffer_A = cols_in_buffer_A + cols_in_buffer_A_now;

            rows_in_buffer_A_now = rows_in_buffer_A_my_initial;
            if(rows_in_buffer_A < rows_in_buffer_A_now)
               rows_in_buffer_A = rows_in_buffer_A_now;

            intNumber = my_pcol/np_rows;
            if(proc_col_min <= my_prow)
               CopyTo = &Buf_to_receive_A[nblk*(na_rows*intNumber - nblk*(intNumber-1)*intNumber/2)];
            else
               CopyTo = &Buf_to_receive_A[nblk*(na_rows*intNumber - nblk*intNumber*(intNumber+1)/2)];
            CopyFrom = Buf_to_send_A;

            Size_receive_A = Size_receive_A + Size_send_A - 1;
         }

         intNumber = ceil((float complex)cols_in_buffer_A_now/(float complex)nblk);
         rows_in_block = rows_in_buffer_A_now;
         for(j = 0; j < intNumber; j++)
         {
            if((j+1)*nblk < cols_in_buffer_A_now)
               cols_in_block = nblk;
            else
               cols_in_block = cols_in_buffer_A_now - j*nblk;

            clacpy("A", &rows_in_block, &cols_in_block, CopyFrom, &rows_in_block, CopyTo, &rows_in_block);

            CopyFrom = CopyFrom + rows_in_block*cols_in_block;
            CopyTo = CopyTo + nblk*(ratio*rows_in_block - nblk*(ratio-1)*ratio/2);
            rows_in_block = rows_in_block - ratio*nblk;
         }
      }
      else
      {
         if(my_prow > 0)
         {
            MPI_Sendrecv(Buf_to_send_A, (int) Size_send_A, MPI_COMPLEX, (int) where_to_send_A, (int) zero, Buf_to_receive_A, (int) Size_U_stored, MPI_COMPLEX, (int) from_where_to_receive_A, (int) zero, row_comm, &status);
            MPI_Get_count(&status, MPI_COMPLEX, &Size_receive_AMPI);
            Size_receive_A = (int) Size_receive_AMPI;

            cols_in_buffer_A = (int)Buf_to_receive_A[Size_receive_A-1];
            if(from_where_to_receive_A <= my_prow)
            {
               rows_in_buffer_A = na_rows;
            }
            else
            {
               rows_in_buffer_A = na_rows - ceil((float complex)(((float complex)from_where_to_receive_A - (float complex)my_prow)/(float complex)np_rows))*nblk;
            }
         }
         else
         {
            Size_receive_A = Size_send_A;
            rows_in_buffer_A = rows_in_buffer_A_my_initial;
            cols_in_buffer_A = cols_in_buffer_A_my_initial;
         }
      }
   }
   if(ratio != 1)
   {
      Buf_to_receive_A[Size_receive_A] = cols_in_buffer_A;
      Buf_to_receive_A[Size_receive_A + 1] = rows_in_buffer_A;
      Size_receive_A = Size_receive_A + 2;
   }
   else
   {
      Buf_to_receive_A[Size_receive_A] = rows_in_buffer_A;
      Size_receive_A = Size_receive_A + 1;
   }

   Size_receive_U = Size_U_skewed;
   U_to_calc = U_stored;

   where_to_send_A = (my_pcol - 1 + np_cols)%np_cols;
   from_where_to_receive_A = (my_pcol + 1)%np_cols;
   where_to_send_U = (my_prow - 1 + np_rows)%np_rows;
   from_where_to_receive_U = (my_prow + 1)%np_rows;
   Curr_pos_in_U_stored = Size_U_skewed;

   for(j = 1; j < np_rows; j++)
   {

      data_ptr = Buf_to_send_A;
      Buf_to_send_A = Buf_to_receive_A;
      Buf_to_receive_A = data_ptr;

      if (j > ToStore)
      {
         data_ptr = Buf_to_send_U;
         Buf_to_send_U = Buf_to_receive_U;
         Buf_to_receive_U = data_ptr;
      }

      Size_send_A = Size_receive_A;
      MPI_Isend(Buf_to_send_A, (int) Size_send_A, MPI_COMPLEX, (int) where_to_send_A, (int) zero, row_comm, &request_A_Send);
      MPI_Irecv(Buf_to_receive_A, (int) (ratio*Size_U_stored), MPI_COMPLEX, (int) from_where_to_receive_A, (int) zero, row_comm, &request_A_Recv);

      Size_send_U = Size_receive_U;
      if (j > ToStore)
      {
         if(j > ToStore + 1)
         {
            MPI_Isend(Buf_to_send_U, (int) Size_send_U, MPI_COMPLEX, (int) where_to_send_U, (int) zero, col_comm, &request_U_Send);
            U_to_calc = Buf_to_send_U;
         }
         else {
     MPI_Isend(U_to_calc, (int) Size_send_U, MPI_COMPLEX, (int) where_to_send_U, (int) zero, col_comm, &request_U_Send);
  }
         MPI_Irecv(Buf_to_receive_U, (int) Size_U_stored, MPI_COMPLEX, (int) from_where_to_receive_U, (int) zero, col_comm, &request_U_Recv);
      }

      rows_in_buffer_U = (int)U_to_calc[Size_receive_U-1];
      row_of_origin_U = (my_pcol + my_prow + np_cols + j - 1)%np_rows;
      if(my_pcol >= row_of_origin_U)
         cols_in_buffer_U = na_cols;
      else
         cols_in_buffer_U = na_cols - nblk;

      cols_in_buffer_A = (int)Buf_to_send_A[Size_receive_A-2];
      rows_in_buffer_A = (int)Buf_to_send_A[Size_receive_A-1];

      col_of_origin_A = np_cols;
      for(i = 0; i < ratio; i++)
      {
         intNumber = (my_pcol + my_prow + i*np_rows + np_cols + j - 1)%np_cols;
         if(intNumber < col_of_origin_A)
            col_of_origin_A = intNumber;
      }

      if (my_pcol >= row_of_origin_U)
         curr_col_loc_res = 0;
      else
         curr_col_loc_res = nblk;

      num_of_blocks_in_U_buffer = ceil((float complex)((float complex)cols_in_buffer_U/(float complex)nblk));
      if(my_pcol >= row_of_origin_U)
         rows_in_block_U = ceil(((float complex)(my_pcol + 1) - (float complex)row_of_origin_U)/(float complex)np_rows)*nblk;
      else
         rows_in_block_U = ratio*nblk;

      U_local_start = U_to_calc;

      for (i = 0; i < num_of_blocks_in_U_buffer; i++)
      {

         curr_col_glob_res = np_cols*nblk*(curr_col_loc_res/nblk) + curr_col_loc_res%nblk + ((np_cols+my_pcol)%np_cols)*nblk;

         Nb = curr_col_glob_res/nblk;
         owner = Nb%np_rows;
         curr_row_loc_res = (Nb/np_rows)*nblk;
         if(my_prow < owner)
            curr_row_loc_res = curr_row_loc_res + nblk;

         curr_row_loc_A = curr_row_loc_res;
         if(col_of_origin_A > my_prow)
            curr_row_loc_A = curr_row_loc_A - nblk;

         rows_in_block = rows_in_buffer_A - curr_row_loc_A;

         curr_col_loc_U = i*nblk;

         if((curr_col_loc_U + nblk) <= cols_in_buffer_U)
            cols_in_block = nblk;
         else
            cols_in_block = cols_in_buffer_U - curr_col_loc_U;

         if(rows_in_block_U > rows_in_buffer_U)
            rows_in_block_U = rows_in_buffer_U;

         A_local_index = curr_row_loc_A;
         A_local_start = &Buf_to_send_A[A_local_index];
         Res_ptr = &Res[curr_col_loc_res*na_rows + curr_row_loc_res];

         LDA_A = rows_in_buffer_A;
         LDA_A_new = LDA_A;
         if ((rows_in_block > 0)&&(cols_in_block > 0))
         {
            U_local_start_curr = U_local_start;

            for (ii = 0; ii < ceil((float complex)rows_in_block_U/(float complex)nblk); ii++)
            {
               if((ii+1)*nblk <= cols_in_buffer_A)
                  rows_in_block_U_curr = nblk;
               else
                  rows_in_block_U_curr = cols_in_buffer_A - ii*nblk;

               if((j == 1)&&(ii == 0)) {
                  cgemm_("N", "N", &rows_in_block, &cols_in_block, &rows_in_block_U_curr, &done, A_local_start, &LDA_A, U_local_start_curr, &rows_in_block_U, &dzero, Res_ptr, &na_rows);
        }
               else {
                  cgemm_("N", "N", &rows_in_block, &cols_in_block, &rows_in_block_U_curr, &done, A_local_start, &LDA_A, U_local_start_curr, &rows_in_block_U, &done, Res_ptr, &na_rows);
               }

               LDA_A_new = LDA_A_new - nblk;

               U_local_start_curr = U_local_start_curr + rows_in_block_U_curr;
               A_local_index = A_local_index - LDA_A + LDA_A*nblk + LDA_A_new;
               A_local_start = &Buf_to_send_A[A_local_index];
               LDA_A = LDA_A_new;
            }
         }

         U_local_start = U_local_start + rows_in_block_U*cols_in_block;
         curr_col_loc_res = curr_col_loc_res + nblk;
         rows_in_block_U = rows_in_block_U + ratio*nblk;
      }

      MPI_Wait(&request_A_Send, &status);
      MPI_Wait(&request_A_Recv, &status);
      MPI_Get_count(&status, MPI_COMPLEX, &Size_receive_AMPI);
      Size_receive_A = (int) Size_receive_AMPI;

      if (j <= ToStore)
      {
         U_to_calc = &U_stored[Curr_pos_in_U_stored];
         Curr_pos_in_U_stored = Curr_pos_in_U_stored + SizesU[j-1];
         Size_receive_U = SizesU[j-1];
      }
      else
      {
         MPI_Wait(&request_U_Send, &status);
         MPI_Wait(&request_U_Recv, &status);
  MPI_Get_count(&status, MPI_COMPLEX, &Size_receive_UMPI);
         Size_receive_U = (int) Size_receive_UMPI;
      }
   }

   if(ToStore < np_rows - 1)
      U_to_calc = Buf_to_receive_U;
   rows_in_buffer_U = (int)U_to_calc[Size_receive_U-1];
   row_of_origin_U = (my_pcol + my_prow + np_cols + j - 1)%np_rows;
   if(my_pcol >= row_of_origin_U)
      cols_in_buffer_U = na_cols;
   else
      cols_in_buffer_U = na_cols - nblk;

   cols_in_buffer_A = (int)Buf_to_receive_A[Size_receive_A-2];
   rows_in_buffer_A = (int)Buf_to_receive_A[Size_receive_A-1];

   col_of_origin_A = np_cols;
   for(i = 0; i < ratio; i++)
   {
      intNumber = (my_pcol + my_prow + i*np_rows + np_cols + np_rows - 1)%np_cols;
      if(intNumber < col_of_origin_A)
         col_of_origin_A = intNumber;
   }

   if (my_pcol >= row_of_origin_U)
      curr_col_loc_res = 0;
   else
      curr_col_loc_res = nblk;

   num_of_blocks_in_U_buffer = ceil((float complex)((float complex)cols_in_buffer_U/(float complex)nblk));
   if(my_pcol >= row_of_origin_U)
      rows_in_block_U = ceil(((float complex)(my_pcol + 1) - (float complex)row_of_origin_U)/(float complex)np_rows)*nblk;
   else
      rows_in_block_U = ratio*nblk;

   U_local_start = U_to_calc;

   for (i = 0; i < num_of_blocks_in_U_buffer; i++)
   {

      curr_col_glob_res = np_cols*nblk*(curr_col_loc_res/nblk) + curr_col_loc_res%nblk + ((np_cols+my_pcol)%np_cols)*nblk;

      Nb = curr_col_glob_res/nblk;
      owner = Nb%np_rows;
      curr_row_loc_res = (Nb/np_rows)*nblk;
      if(my_prow < owner)
         curr_row_loc_res = curr_row_loc_res + nblk;

      curr_row_loc_A = curr_row_loc_res;
      if(col_of_origin_A > my_prow)
         curr_row_loc_A = curr_row_loc_A - nblk;

      rows_in_block = rows_in_buffer_A - curr_row_loc_A;

      curr_col_loc_U = i*nblk;

      if((curr_col_loc_U + nblk) <= cols_in_buffer_U)
         cols_in_block = nblk;
      else
         cols_in_block = cols_in_buffer_U - curr_col_loc_U;

      if(rows_in_block_U > rows_in_buffer_U)
         rows_in_block_U = rows_in_buffer_U;

      A_local_index = curr_row_loc_A;
      A_local_start = &Buf_to_receive_A[A_local_index];
      Res_ptr = &Res[curr_col_loc_res*na_rows + curr_row_loc_res];
      LDA_A = rows_in_buffer_A;
      LDA_A_new = LDA_A;
      if ((rows_in_block > 0) &&(cols_in_block > 0))
      {
         U_local_start_curr = U_local_start;

         for (ii = 0; ii < ceil((float complex)rows_in_block_U/(float complex)nblk); ii++)
         {
            if((ii+1)*nblk <= cols_in_buffer_A)
               rows_in_block_U_curr = nblk;
            else
               rows_in_block_U_curr = cols_in_buffer_A - ii*nblk;

            if((j == 1)&&(ii == 0)) {
               cgemm_("N", "N", &rows_in_block, &cols_in_block, &rows_in_block_U_curr, &done, A_local_start, &LDA_A, U_local_start_curr, &rows_in_block_U, &dzero, Res_ptr, &na_rows);
     }
            else {
               cgemm_("N", "N", &rows_in_block, &cols_in_block, &rows_in_block_U_curr, &done, A_local_start, &LDA_A, U_local_start_curr, &rows_in_block_U, &done, Res_ptr, &na_rows);
     }

            LDA_A_new = LDA_A_new - nblk;

            U_local_start_curr = U_local_start_curr + rows_in_block_U_curr;
            A_local_index = A_local_index - (LDA_A - rows_in_block) + LDA_A*nblk + LDA_A_new - rows_in_block;
            A_local_start = &Buf_to_receive_A[A_local_index];
            LDA_A = LDA_A_new;
         }
      }

      U_local_start = U_local_start + rows_in_block_U*cols_in_block;
      curr_col_loc_res = curr_col_loc_res + nblk;
      rows_in_block_U = rows_in_block_U + ratio*nblk;
   }

   pctranc(&na, &na, &done, Res, &one, &one, a_desc, &dzero, M, &one, &one, a_desc);
   pclacpy("U", &na, &na, M, &one, &one, a_desc, Res, &one, &one, a_desc);

   free(Buf_to_send_A);
   free(Buf_to_receive_A);
   free(Buf_to_send_U);
   free(Buf_to_receive_U);
   free(M);
   free(M_T);
   if(ratio != 1)
      free(Buf_A);
   free(U_stored);
   free(SizesU);
}

void cannons_reduction_c_fc(float complex* A, float complex* U, int local_rowsCast, int local_colsCast,
                         int* a_desc, float complex *Res, int ToStore, int row_comm, int col_comm)
{
  int local_rows, local_cols;
  local_rows = (int) local_rowsCast;
  local_cols = (int) local_colsCast;

  MPI_Comm c_row_comm = MPI_Comm_f2c(row_comm);
  MPI_Comm c_col_comm = MPI_Comm_f2c(col_comm);

  int my_prow, my_pcol, np_rows, np_cols;
  int my_prowMPI, my_pcolMPI, np_rowsMPI, np_colsMPI;

  MPI_Comm_rank(c_row_comm, &my_prowMPI);
  MPI_Comm_size(c_row_comm, &np_rowsMPI);
  MPI_Comm_rank(c_col_comm, &my_pcolMPI);
  MPI_Comm_size(c_col_comm, &np_colsMPI);

  my_prow = (int) my_prowMPI;
  my_pcol = (int) my_pcolMPI;
  np_rows = (int) np_rowsMPI;
  np_cols = (int) np_colsMPI;

  cannons_reduction_fc(A, U, np_rows, np_cols, my_prow, my_pcol, a_desc, Res, ToStore, c_col_comm, c_row_comm);
}
void cannons_triang_rectangular_fc(float complex* U, float complex* B, int np_rows, int np_cols, int my_prow, int my_pcol, int* U_desc, int* b_desc, float complex *Res, MPI_Comm row_comm, MPI_Comm col_comm)
{
   int na, nb, nblk, width, na_rows, na_cols, nb_cols, cols_in_buffer_U_my_initial, cols_in_buffer_U, rows_in_buffer_U, Size_receive_U_now, rows_in_buffer_U_now, cols_in_buffer_U_now, rows_in_buffer_U_my_initial;

   int Size_receive_U_nowMPI, Size_receive_UMPI, Size_receive_BMPI;
   int i, j, Size_send_U, Size_receive_U, Size_send_B, Size_receive_B, intNumber, Buf_rows, Buf_cols_U, Buf_cols_B, curr_rows, num_of_iters, cols_in_buffer, rows_in_block, curr_col_loc, cols_in_block, num_of_blocks_in_U_buffer, col_of_origin_U, b_rows_mult, b_cols_mult;

   float complex *Buf_to_send_U, *Buf_to_receive_U, *Buf_to_send_B, *Buf_to_receive_B, *Buf_U, *PosBuff;

   int where_to_send_U, from_where_to_receive_U, where_to_send_B, from_where_to_receive_B, last_proc_col_B, last_proc_row_B, n, Size_U_stored, proc_col_min;

   float complex *U_local_start, *Buf_pos, *B_local_start, *double_ptr, *CopyTo, *CopyFrom;

   int ratio;

   MPI_Status status;

   int one = 1;
   int zero = 0;
   float complex done = 1.0;
   float complex dzero = 0.0;

   na = U_desc[2];
   nblk = U_desc[4];
   nb = b_desc[3];

   na_rows = numroc(&na, &nblk, &my_prow, &zero, &np_rows);
   na_cols = numroc(&na, &nblk, &my_pcol, &zero, &np_cols);
   nb_cols = numroc(&nb, &nblk, &my_pcol, &zero, &np_cols);

   MPI_Request request_U_Recv;
   MPI_Request request_U_Send;
   MPI_Request request_B_Recv;
   MPI_Request request_B_Send;

   last_proc_col_B = ((nb-1)/nblk) % np_cols;
   last_proc_row_B = ((na-1)/nblk) % np_rows;

    if(nb%nblk == 0)
      if(my_pcol <= last_proc_col_B)
         Buf_cols_B = nb_cols;
      else
         Buf_cols_B = nb_cols + nblk;
   else
      if(my_pcol < last_proc_col_B)
         Buf_cols_B = nb_cols;
      else if(my_pcol > last_proc_col_B)
         Buf_cols_B = nb_cols + nblk;
      else
         Buf_cols_B = nb_cols + nblk - nb_cols%nblk;

   if(na%nblk == 0)
      if(my_prow <= last_proc_row_B)
         Buf_rows = na_rows;
      else
         Buf_rows = na_rows + nblk;
   else
      if(my_prow < last_proc_row_B)
         Buf_rows = na_rows;
      else if(my_prow > last_proc_row_B)
         Buf_rows = na_rows + nblk;
      else
         Buf_rows = na_rows + nblk - na_rows%nblk;

   ratio = np_cols/np_rows;

   intNumber = ceil((float complex)na/(float complex)(np_cols*nblk));
   Size_U_stored = ratio*nblk*nblk*intNumber*(intNumber+1)/2 + 2;

   Buf_to_send_U = malloc(ratio*Size_U_stored*sizeof(float complex));
   Buf_to_receive_U = malloc(ratio*Size_U_stored*sizeof(float complex));
   Buf_to_send_B = malloc(Buf_cols_B*Buf_rows*sizeof(float complex));
   Buf_to_receive_B = malloc(Buf_cols_B*Buf_rows*sizeof(float complex));
   if(ratio != 1)
      Buf_U = malloc(Size_U_stored*sizeof(float complex));

   for(i = 0; i < na_rows*nb_cols; i++)
     Res[i] = 0;

   if((ratio != 1)||(my_prow != 0))
      Buf_pos = Buf_to_send_U;
   else
      Buf_pos = Buf_to_receive_U;

   if(my_pcol >= my_prow)
      curr_col_loc = 0;
   else
      curr_col_loc = 1;

   num_of_iters = ceil((float complex)na_cols/(float complex)nblk);
   num_of_iters = num_of_iters - curr_col_loc;
   curr_col_loc = curr_col_loc*nblk;

   if(my_pcol >= my_prow )
      rows_in_block = ceil(((float complex)(my_pcol + 1) - (float complex)my_prow)/(float complex)np_rows)*nblk;
   else
      rows_in_block = ratio*nblk;
   cols_in_buffer_U_my_initial = 0;
   Size_send_U = 0;
   for(i = 0; i < num_of_iters; i++)
   {
      if(rows_in_block > na_rows)
         rows_in_block = na_rows;

      if ((na_cols - curr_col_loc) < nblk)
         cols_in_block = na_cols - curr_col_loc;
      else
         cols_in_block = nblk;

      if((rows_in_block > 0)&&(cols_in_block > 0))
      {
         double_ptr = &U[curr_col_loc*na_rows];
         clacpy("A", &rows_in_block, &cols_in_block, double_ptr, &na_rows, Buf_pos, &rows_in_block);
         Buf_pos = Buf_pos + rows_in_block*cols_in_block;
         Size_send_U = Size_send_U + rows_in_block*cols_in_block;
         cols_in_buffer_U_my_initial = cols_in_buffer_U_my_initial + cols_in_block;
      }
      curr_col_loc = curr_col_loc + nblk;
      rows_in_block = rows_in_block + ratio*nblk;
   }
   rows_in_buffer_U_my_initial = rows_in_block - ratio*nblk;
   *Buf_pos = (float complex)cols_in_buffer_U_my_initial;
   Buf_pos = Buf_pos + 1;
   *Buf_pos = (float complex)rows_in_buffer_U_my_initial;
   Size_send_U = Size_send_U + 2;

   proc_col_min = np_cols;
   for(i = 0; i < ratio; i++)
   {
      from_where_to_receive_U = (my_pcol + my_prow + i*np_rows)%np_cols;
      if(from_where_to_receive_U < proc_col_min)
         proc_col_min = from_where_to_receive_U;
   }

   Size_receive_U = 0;
   cols_in_buffer_U = 0;
   rows_in_buffer_U = 0;
   for(i = 0; i < ratio; i++)
   {
      where_to_send_U = (my_pcol - my_prow - i*np_rows + np_cols)%np_cols;
      from_where_to_receive_U = (my_pcol + my_prow + i*np_rows)%np_cols;

      if(ratio != 1)
      {
         if(where_to_send_U != my_pcol)
         {
            MPI_Sendrecv(Buf_to_send_U, (int) Size_send_U, MPI_COMPLEX, (int) where_to_send_U, 0, Buf_U, (int) Size_U_stored, MPI_COMPLEX, (int) from_where_to_receive_U, 0, row_comm, &status);
            MPI_Get_count(&status, MPI_COMPLEX, &Size_receive_U_nowMPI);
            Size_receive_U_now = (int) Size_receive_U_nowMPI;
            Size_receive_U = Size_receive_U + Size_receive_U_now - 2;

            cols_in_buffer_U_now = Buf_U[Size_receive_U_now - 2];
            cols_in_buffer_U = cols_in_buffer_U + cols_in_buffer_U_now;
            rows_in_buffer_U_now = Buf_U[Size_receive_U_now - 1];

            if(rows_in_buffer_U < rows_in_buffer_U_now)
               rows_in_buffer_U = rows_in_buffer_U_now;

            intNumber = from_where_to_receive_U/np_rows;
            if(proc_col_min >= my_prow)
               CopyTo = &Buf_to_receive_U[nblk*nblk*intNumber*(intNumber + 1)/2];
            else
               if(from_where_to_receive_U < my_prow)
                  CopyTo = &Buf_to_receive_U[nblk*nblk*ratio*(ratio - 1)/2];
               else
                  CopyTo = &Buf_to_receive_U[nblk*nblk*intNumber*(intNumber - 1)/2];
            CopyFrom = Buf_U;
         }
         else
         {
            cols_in_buffer_U_now = cols_in_buffer_U_my_initial;
            cols_in_buffer_U = cols_in_buffer_U + cols_in_buffer_U_now;

            rows_in_buffer_U_now = rows_in_buffer_U_my_initial;
            if(rows_in_buffer_U < rows_in_buffer_U_now)
               rows_in_buffer_U = rows_in_buffer_U_now;

            intNumber = my_pcol/np_rows;
            if(proc_col_min >= my_prow)
               CopyTo = &Buf_to_receive_U[nblk*nblk*intNumber*(intNumber + 1)/2];
            else
               if(my_pcol < my_prow)
                  CopyTo = &Buf_to_receive_U[nblk*nblk*ratio*(ratio - 1)/2];
               else
                  CopyTo = &Buf_to_receive_U[nblk*nblk*intNumber*(intNumber - 1)/2];
            CopyFrom = Buf_to_send_U;
            Size_receive_U = Size_receive_U + Size_send_U - 2;
         }

         intNumber = ceil((float complex)cols_in_buffer_U_now/(float complex)nblk);
         if(from_where_to_receive_U >= my_prow)
            rows_in_block = ceil(((float complex)(from_where_to_receive_U + 1) - (float complex)my_prow)/(float complex)np_rows)*nblk;
         else
            rows_in_block = ratio*nblk;
         for(j = 0; j < intNumber; j++)
         {
            if((j+1)*nblk < cols_in_buffer_U_now)
               cols_in_block = nblk;
            else
               cols_in_block = cols_in_buffer_U_now - j*nblk;

            clacpy("A", &rows_in_block, &cols_in_block, CopyFrom, &rows_in_block, CopyTo, &rows_in_block);

            CopyFrom = CopyFrom + rows_in_block*cols_in_block;
            CopyTo = CopyTo + ratio*rows_in_block*nblk + nblk*nblk*ratio*(ratio-1)/2;
            rows_in_block = rows_in_block + ratio*nblk;
            if(rows_in_block > rows_in_buffer_U_now)
               rows_in_block = rows_in_buffer_U_now;
         }
      }
      else
      {
         if(my_prow > 0)
         {
            MPI_Sendrecv(Buf_to_send_U, (int) Size_send_U, MPI_COMPLEX, (int) where_to_send_U, 0, Buf_to_receive_U, (int) Size_U_stored, MPI_COMPLEX, (int) from_where_to_receive_U, 0, row_comm, &status);
            MPI_Get_count(&status, MPI_COMPLEX, &Size_receive_UMPI);
            Size_receive_U = (int) Size_receive_UMPI;

            cols_in_buffer_U = (int)Buf_to_receive_U[Size_receive_U-2];
            rows_in_buffer_U = (int)Buf_to_receive_U[Size_receive_U-1];
         }
         else
         {
            Size_receive_U = Size_send_U;
            rows_in_buffer_U = rows_in_buffer_U_my_initial;
            cols_in_buffer_U = cols_in_buffer_U_my_initial;
         }
      }
   }
   if(ratio != 1)
   {
      Buf_to_receive_U[Size_receive_U] = cols_in_buffer_U;
      Buf_to_receive_U[Size_receive_U + 1] = rows_in_buffer_U;
      Size_receive_U = Size_receive_U + 2;
   }

   if(my_pcol > 0)
   {
      where_to_send_B = (my_prow - my_pcol + np_cols)%np_rows;
      from_where_to_receive_B = (my_pcol + my_prow)%np_rows;

      if(where_to_send_B != my_prow)
      {

         clacpy("A", &na_rows, &nb_cols, B, &na_rows, Buf_to_send_B, &na_rows);
         MPI_Sendrecv(Buf_to_send_B, (int) (nb_cols*na_rows), MPI_COMPLEX, (int) where_to_send_B, 0, Buf_to_receive_B, (int) (nb_cols*Buf_rows), MPI_COMPLEX, (int) from_where_to_receive_B, 0, col_comm, &status);
         MPI_Get_count(&status, MPI_COMPLEX, &Size_receive_BMPI);
         Size_receive_B = (int) Size_receive_BMPI;
         Size_receive_B = Size_receive_B/nb_cols;

      }
      else
      {
         clacpy("A", &na_rows, &nb_cols, B, &na_rows, Buf_to_receive_B, &na_rows);
         Size_receive_B = na_rows;
      }
   }
   else
   {
      clacpy("A", &na_rows, &nb_cols, B, &na_rows, Buf_to_receive_B, &na_rows);
      Size_receive_B = na_rows;
   }

   where_to_send_U = (my_pcol - 1 + np_cols)%np_cols;
   from_where_to_receive_U = (my_pcol + 1)%np_cols;
   where_to_send_B = (my_prow - 1 + np_rows)%np_rows;
   from_where_to_receive_B = (my_prow + 1)%np_rows;

   for(i = 1; i < np_rows; i++)
   {

      double_ptr = Buf_to_send_U;
      Buf_to_send_U = Buf_to_receive_U;
      Buf_to_receive_U = double_ptr;

      double_ptr = Buf_to_send_B;
      Buf_to_send_B = Buf_to_receive_B;
      Buf_to_receive_B = double_ptr;

      Size_send_U = Size_receive_U;
      Size_send_B = Size_receive_B;

      MPI_Isend(Buf_to_send_U, (int) Size_send_U, MPI_COMPLEX, (int) where_to_send_U, 0, row_comm, &request_U_Send);
      MPI_Irecv(Buf_to_receive_U, (int) (ratio*Size_U_stored), MPI_COMPLEX, (int) from_where_to_receive_U, 0, row_comm, &request_U_Recv);

      MPI_Isend(Buf_to_send_B, (int) (Size_send_B*nb_cols), MPI_COMPLEX, (int) where_to_send_B, 0, col_comm, &request_B_Send);
      MPI_Irecv(Buf_to_receive_B, (int) (Buf_rows*nb_cols), MPI_COMPLEX, (int) from_where_to_receive_B, 0, col_comm, &request_B_Recv);

      cols_in_buffer_U = (int)Buf_to_send_U[Size_receive_U-2];
      rows_in_buffer_U = (int)Buf_to_send_U[Size_receive_U-1];

      proc_col_min = np_cols;
      for(j = 0; j < ratio; j++)
      {
         col_of_origin_U = (my_pcol + my_prow + i - 1 + j*np_rows)%np_cols;
         if(col_of_origin_U < proc_col_min)
            proc_col_min = col_of_origin_U;
      }
      col_of_origin_U = proc_col_min;

      num_of_blocks_in_U_buffer = ceil((float complex)cols_in_buffer_U/(float complex)nblk);

      if (col_of_origin_U >= my_prow)
         B_local_start = Buf_to_send_B;
      else
         B_local_start = Buf_to_send_B + nblk;

      U_local_start = Buf_to_send_U;

      for(j = 0; j < num_of_blocks_in_U_buffer; j++)
      {
         curr_rows = (j+1)*nblk;
         if (curr_rows > rows_in_buffer_U)
            curr_rows = rows_in_buffer_U;

         if((j+1)*nblk <= cols_in_buffer_U)
            b_rows_mult = nblk;
         else
            b_rows_mult = cols_in_buffer_U - j*nblk;

         cgemm_("N", "N", &curr_rows, &nb_cols, &b_rows_mult, &done, U_local_start, &curr_rows, B_local_start, &Size_receive_B, &done, Res, &na_rows);

         U_local_start = U_local_start + nblk*curr_rows;
         B_local_start = B_local_start + nblk;
      }

      MPI_Wait(&request_U_Send, &status);
      MPI_Wait(&request_U_Recv, &status);
      MPI_Get_count(&status, MPI_COMPLEX, &Size_receive_UMPI);
      Size_receive_U = (int) Size_receive_UMPI;

      MPI_Wait(&request_B_Send, &status);
      MPI_Wait(&request_B_Recv, &status);
      MPI_Get_count(&status, MPI_COMPLEX, &Size_receive_BMPI);
      Size_receive_B = (int) Size_receive_BMPI;
      Size_receive_B = (int) Size_receive_B / nb_cols;

   }

   cols_in_buffer_U = (int)Buf_to_receive_U[Size_receive_U-2];
   rows_in_buffer_U = (int)Buf_to_receive_U[Size_receive_U-1];

   proc_col_min = np_cols;
   for(j = 0; j < ratio; j++)
   {
      col_of_origin_U = (my_pcol + my_prow + np_rows - 1 + j*np_rows)%np_cols;
      if(col_of_origin_U < proc_col_min)
         proc_col_min = col_of_origin_U;
   }
   col_of_origin_U = proc_col_min;

   num_of_blocks_in_U_buffer = ceil((float complex)cols_in_buffer_U/(float complex)nblk);

   if (col_of_origin_U >= my_prow)
      B_local_start = Buf_to_receive_B;
   else
      B_local_start = Buf_to_receive_B + nblk;

   U_local_start = Buf_to_receive_U;

   for(j = 0; j < num_of_blocks_in_U_buffer; j++)
   {
      curr_rows = (j+1)*nblk;
      if (curr_rows > rows_in_buffer_U)
         curr_rows = rows_in_buffer_U;

      if((j+1)*nblk <= cols_in_buffer_U)
         b_rows_mult = nblk;
      else
         b_rows_mult = cols_in_buffer_U - j*nblk;

      cgemm_("N", "N", &curr_rows, &nb_cols, &b_rows_mult, &done, U_local_start, &curr_rows, B_local_start, &Size_receive_B, &done, Res, &na_rows);

      U_local_start = U_local_start + nblk*curr_rows;
      B_local_start = B_local_start + nblk;
   }

   free(Buf_to_send_U);
   free(Buf_to_receive_U);
   free(Buf_to_send_B);
   free(Buf_to_receive_B);
   if(ratio != 1)
      free(Buf_U);
}

void cannons_triang_rectangular_c_fc(float complex* U, float complex* B, int local_rowsCast, int local_colsCast,
                                    int* u_desc, int* b_desc, float complex *Res, int row_comm, int col_comm)
{
  int local_rows, local_cols;

  local_rows = (int) local_rowsCast;
  local_cols = (int) local_colsCast;

  MPI_Comm c_row_comm = MPI_Comm_f2c(row_comm);
  MPI_Comm c_col_comm = MPI_Comm_f2c(col_comm);

  int my_prow, my_pcol, np_rows, np_cols;
  int my_prowMPI, my_pcolMPI, np_rowsMPI, np_colsMPI;

  MPI_Comm_rank(c_row_comm, &my_prowMPI);
  MPI_Comm_size(c_row_comm, &np_rowsMPI);
  MPI_Comm_rank(c_col_comm, &my_pcolMPI);
  MPI_Comm_size(c_col_comm, &np_colsMPI);

  my_prow = (int) my_prowMPI;
  my_pcol = (int) my_pcolMPI;
  np_rows = (int) np_rowsMPI;
  np_cols = (int) np_colsMPI;

  cannons_triang_rectangular_fc(U, B, np_rows, np_cols, my_prow, my_pcol, u_desc, b_desc, Res, c_col_comm, c_row_comm);
}
void cannons_reduction_c_fc(float complex* A, float complex* U, int local_rowsCast, int local_colsCast, int* a_desc,
                         float complex *Res, int ToStore, int row_comm, int col_comm);
void cannons_triang_rectangular_c_fc(float complex* U, float complex* B, int local_rowsCast, int local_colsCast,
                                    int* u_desc, int* b_desc, float complex *Res, int row_comm, int col_comm);
