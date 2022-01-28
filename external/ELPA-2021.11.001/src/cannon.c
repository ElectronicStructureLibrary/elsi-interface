














 








 

 


 

 









































































 



 

















 











































































 


 


 

 







 


 

 


 

 

 

 




 


 



















 






 


 


 

 
















 


 



 

 



 





 



 


 


 


 






 


 


 


 


 
 










 







 



 



 



 

 



 

 






 




 


 


 

 


 









 



 




 

 



 


 


 




 


 





 


 





 






































 
typedef long ptrdiff_t;

typedef unsigned long size_t;


 

 
typedef int wchar_t;



typedef struct {
  long long __clang_max_align_nonce1
      __attribute__((__aligned__(__alignof__(long long))));
  long double __clang_max_align_nonce2
      __attribute__((__aligned__(__alignof__(long double))));
} max_align_t;



 

















 



 


 


 

 
typedef unsigned char __u_char;
typedef unsigned short int __u_short;
typedef unsigned int __u_int;
typedef unsigned long int __u_long;

 
typedef signed char __int8_t;
typedef unsigned char __uint8_t;
typedef signed short int __int16_t;
typedef unsigned short int __uint16_t;
typedef signed int __int32_t;
typedef unsigned int __uint32_t;
typedef signed long int __int64_t;
typedef unsigned long int __uint64_t;

 
typedef long int __quad_t;
typedef unsigned long int __u_quad_t;





























 

 
















 




 

 




 

 

 




typedef unsigned long int __dev_t;	 
typedef unsigned int __uid_t;	 
typedef unsigned int __gid_t;	 
typedef unsigned long int __ino_t;	 
typedef unsigned long int __ino64_t;	 
typedef unsigned int __mode_t;	 
typedef unsigned long int __nlink_t;	 
typedef long int __off_t;	 
typedef long int __off64_t;	 
typedef int __pid_t;	 
typedef struct { int __val[2]; } __fsid_t;	 
typedef long int __clock_t;	 
typedef unsigned long int __rlim_t;	 
typedef unsigned long int __rlim64_t;	 
typedef unsigned int __id_t;		 
typedef long int __time_t;	 
typedef unsigned int __useconds_t;  
typedef long int __suseconds_t;  

typedef int __daddr_t;	 
typedef int __key_t;	 

 
typedef int __clockid_t;

 
typedef void * __timer_t;

 
typedef long int __blksize_t;

 

 
typedef long int __blkcnt_t;
typedef long int __blkcnt64_t;

 
typedef unsigned long int __fsblkcnt_t;
typedef unsigned long int __fsblkcnt64_t;

 
typedef unsigned long int __fsfilcnt_t;
typedef unsigned long int __fsfilcnt64_t;

 
typedef long int __fsword_t;

typedef long int __ssize_t;  

 
typedef long int __syscall_slong_t;
 
typedef unsigned long int __syscall_ulong_t;


 
typedef __off64_t __loff_t;	 
typedef __quad_t *__qaddr_t;
typedef char *__caddr_t;

 
typedef long int __intptr_t;

 
typedef unsigned int __socklen_t;






 
struct _IO_FILE;


 
typedef struct _IO_FILE FILE;






 
typedef struct _IO_FILE __FILE;






























 



 


 
































 



 

 






 















 




 




 
typedef struct
{
  int __count;
  union
  {
    unsigned int __wch;
    char __wchb[4];
  } __value;		 
} __mbstate_t;



 



 
typedef struct
{
  __off_t __pos;
  __mbstate_t __state;
} _G_fpos_t;
typedef struct
{
  __off64_t __pos;
  __mbstate_t __state;
} _G_fpos64_t;


 



 


 

 








 























 


typedef __builtin_va_list va_list;



 


 
typedef __builtin_va_list __gnuc_va_list;










 



 


struct _IO_jump_t;  struct _IO_FILE;

 
typedef void _IO_lock_t;


 

struct _IO_marker {
  struct _IO_marker *_next;
  struct _IO_FILE *_sbuf;
  
 
   
  int _pos;
};

 
enum __codecvt_result
{
  __codecvt_ok,
  __codecvt_partial,
  __codecvt_error,
  __codecvt_noconv
};


struct _IO_FILE {
  int _flags;		 

   
   
  char* _IO_read_ptr;	 
  char* _IO_read_end;	 
  char* _IO_read_base;	 
  char* _IO_write_base;	 
  char* _IO_write_ptr;	 
  char* _IO_write_end;	 
  char* _IO_buf_base;	 
  char* _IO_buf_end;	 
   
  char *_IO_save_base;  
  char *_IO_backup_base;   
  char *_IO_save_end;  

  struct _IO_marker *_markers;

  struct _IO_FILE *_chain;

  int _fileno;
  int _flags2;
  __off_t _old_offset;  

   
  unsigned short _cur_column;
  signed char _vtable_offset;
  char _shortbuf[1];

   

  _IO_lock_t *_lock;
  __off64_t _offset;
  void *__pad1;
  void *__pad2;
  void *__pad3;
  void *__pad4;
  size_t __pad5;
  int _mode;
   
  char _unused2[15 * sizeof (int) - 4 * sizeof (void *) - sizeof (size_t)];
};

typedef struct _IO_FILE _IO_FILE;

struct _IO_FILE_plus;

extern struct _IO_FILE_plus _IO_2_1_stdin_;
extern struct _IO_FILE_plus _IO_2_1_stdout_;
extern struct _IO_FILE_plus _IO_2_1_stderr_;


 


 
typedef __ssize_t __io_read_fn (void *__cookie, char *__buf, size_t __nbytes);






 
typedef __ssize_t __io_write_fn (void *__cookie, const char *__buf,
				 size_t __n);






 
typedef int __io_seek_fn (void *__cookie, __off64_t *__pos, int __w);

 
typedef int __io_close_fn (void *__cookie);





extern int __underflow (_IO_FILE *);
extern int __uflow (_IO_FILE *);
extern int __overflow (_IO_FILE *, int);





extern int _IO_getc (_IO_FILE *__fp);
extern int _IO_putc (int __c, _IO_FILE *__fp);
extern int _IO_feof (_IO_FILE *__fp) __attribute__ ((__nothrow__ , __leaf__));
extern int _IO_ferror (_IO_FILE *__fp) __attribute__ ((__nothrow__ , __leaf__));

extern int _IO_peekc_locked (_IO_FILE *__fp);

 

extern void _IO_flockfile (_IO_FILE *) __attribute__ ((__nothrow__ , __leaf__));
extern void _IO_funlockfile (_IO_FILE *) __attribute__ ((__nothrow__ , __leaf__));
extern int _IO_ftrylockfile (_IO_FILE *) __attribute__ ((__nothrow__ , __leaf__));


extern int _IO_vfscanf (_IO_FILE * __restrict, const char * __restrict,
			__gnuc_va_list, int *__restrict);
extern int _IO_vfprintf (_IO_FILE *__restrict, const char *__restrict,
			 __gnuc_va_list);
extern __ssize_t _IO_padn (_IO_FILE *, int, __ssize_t);
extern size_t _IO_sgetn (_IO_FILE *, void *, size_t);

extern __off64_t _IO_seekoff (_IO_FILE *, __off64_t, int, int);
extern __off64_t _IO_seekpos (_IO_FILE *, __off64_t, int);

extern void _IO_free_backup_area (_IO_FILE *) __attribute__ ((__nothrow__ , __leaf__));





typedef __gnuc_va_list va_list;

typedef __off_t off_t;

typedef __ssize_t ssize_t;

 

typedef _G_fpos_t fpos_t;


 


 



 



 


 










 















 







 
extern struct _IO_FILE *stdin;		 
extern struct _IO_FILE *stdout;		 
extern struct _IO_FILE *stderr;		 
 


 
extern int remove (const char *__filename) __attribute__ ((__nothrow__ , __leaf__));
 
extern int rename (const char *__old, const char *__new) __attribute__ ((__nothrow__ , __leaf__));


 
extern int renameat (int __oldfd, const char *__old, int __newfd,
		     const char *__new) __attribute__ ((__nothrow__ , __leaf__));





 
extern FILE *tmpfile (void) ;


 
extern char *tmpnam (char *__s) __attribute__ ((__nothrow__ , __leaf__)) ;



 
extern char *tmpnam_r (char *__s) __attribute__ ((__nothrow__ , __leaf__)) ;








 
extern char *tempnam (const char *__dir, const char *__pfx)
     __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__malloc__)) ;






 
extern int fclose (FILE *__stream);



 
extern int fflush (FILE *__stream);







 
extern int fflush_unlocked (FILE *__stream);







 
extern FILE *fopen (const char *__restrict __filename,
		    const char *__restrict __modes) ;



 
extern FILE *freopen (const char *__restrict __filename,
		      const char *__restrict __modes,
		      FILE *__restrict __stream) ;


 
extern FILE *fdopen (int __fd, const char *__modes) __attribute__ ((__nothrow__ , __leaf__)) ;


 
extern FILE *fmemopen (void *__s, size_t __len, const char *__modes)
  __attribute__ ((__nothrow__ , __leaf__)) ;



 
extern FILE *open_memstream (char **__bufloc, size_t *__sizeloc) __attribute__ ((__nothrow__ , __leaf__)) ;




 
extern void setbuf (FILE *__restrict __stream, char *__restrict __buf) __attribute__ ((__nothrow__ , __leaf__));


 
extern int setvbuf (FILE *__restrict __stream, char *__restrict __buf,
		    int __modes, size_t __n) __attribute__ ((__nothrow__ , __leaf__));



 
extern void setbuffer (FILE *__restrict __stream, char *__restrict __buf,
		       size_t __size) __attribute__ ((__nothrow__ , __leaf__));

 
extern void setlinebuf (FILE *__stream) __attribute__ ((__nothrow__ , __leaf__));






 
extern int fprintf (FILE *__restrict __stream,
		    const char *__restrict __format, ...);



 
extern int printf (const char *__restrict __format, ...);
 
extern int sprintf (char *__restrict __s,
		    const char *__restrict __format, ...) __attribute__ ((__nothrow__));




 
extern int vfprintf (FILE *__restrict __s, const char *__restrict __format,
		     __gnuc_va_list __arg);



 
extern int vprintf (const char *__restrict __format, __gnuc_va_list __arg);
 
extern int vsprintf (char *__restrict __s, const char *__restrict __format,
		     __gnuc_va_list __arg) __attribute__ ((__nothrow__));



 
extern int snprintf (char *__restrict __s, size_t __maxlen,
		     const char *__restrict __format, ...)
     __attribute__ ((__nothrow__)) __attribute__ ((__format__ (__printf__, 3, 4)));

extern int vsnprintf (char *__restrict __s, size_t __maxlen,
		      const char *__restrict __format, __gnuc_va_list __arg)
     __attribute__ ((__nothrow__)) __attribute__ ((__format__ (__printf__, 3, 0)));



 
extern int vdprintf (int __fd, const char *__restrict __fmt,
		     __gnuc_va_list __arg)
     __attribute__ ((__format__ (__printf__, 2, 0)));
extern int dprintf (int __fd, const char *__restrict __fmt, ...)
     __attribute__ ((__format__ (__printf__, 2, 3)));






 
extern int fscanf (FILE *__restrict __stream,
		   const char *__restrict __format, ...) ;



 
extern int scanf (const char *__restrict __format, ...) ;
 
extern int sscanf (const char *__restrict __s,
		   const char *__restrict __format, ...) __attribute__ ((__nothrow__ , __leaf__));



 
extern int fscanf (FILE *__restrict __stream, const char *__restrict __format, ...) __asm__ ("" "__isoc99_fscanf") ;
extern int scanf (const char *__restrict __format, ...) __asm__ ("" "__isoc99_scanf") ;
extern int sscanf (const char *__restrict __s, const char *__restrict __format, ...) __asm__ ("" "__isoc99_sscanf") __attribute__ ((__nothrow__ , __leaf__));







 
extern int vfscanf (FILE *__restrict __s, const char *__restrict __format,
		    __gnuc_va_list __arg)
     __attribute__ ((__format__ (__scanf__, 2, 0))) ;




 
extern int vscanf (const char *__restrict __format, __gnuc_va_list __arg)
     __attribute__ ((__format__ (__scanf__, 1, 0))) ;

 
extern int vsscanf (const char *__restrict __s,
		    const char *__restrict __format, __gnuc_va_list __arg)
     __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__format__ (__scanf__, 2, 0)));



 
extern int vfscanf (FILE *__restrict __s, const char *__restrict __format, __gnuc_va_list __arg) __asm__ ("" "__isoc99_vfscanf")
     __attribute__ ((__format__ (__scanf__, 2, 0))) ;
extern int vscanf (const char *__restrict __format, __gnuc_va_list __arg) __asm__ ("" "__isoc99_vscanf")
     __attribute__ ((__format__ (__scanf__, 1, 0))) ;
extern int vsscanf (const char *__restrict __s, const char *__restrict __format, __gnuc_va_list __arg) __asm__ ("" "__isoc99_vsscanf") __attribute__ ((__nothrow__ , __leaf__))
     __attribute__ ((__format__ (__scanf__, 2, 0)));








 
extern int fgetc (FILE *__stream);
extern int getc (FILE *__stream);




 
extern int getchar (void);



 




 
extern int getc_unlocked (FILE *__stream);
extern int getchar_unlocked (void);






 
extern int fgetc_unlocked (FILE *__stream);









 
extern int fputc (int __c, FILE *__stream);
extern int putc (int __c, FILE *__stream);




 
extern int putchar (int __c);



 






 
extern int fputc_unlocked (int __c, FILE *__stream);




 
extern int putc_unlocked (int __c, FILE *__stream);
extern int putchar_unlocked (int __c);


 
extern int getw (FILE *__stream);

 
extern int putw (int __w, FILE *__stream);






 
extern char *fgets (char *__restrict __s, int __n, FILE *__restrict __stream)
     ;














 
extern __ssize_t __getdelim (char **__restrict __lineptr,
			       size_t *__restrict __n, int __delimiter,
			       FILE *__restrict __stream) ;
extern __ssize_t getdelim (char **__restrict __lineptr,
			     size_t *__restrict __n, int __delimiter,
			     FILE *__restrict __stream) ;






 
extern __ssize_t getline (char **__restrict __lineptr,
			    size_t *__restrict __n,
			    FILE *__restrict __stream) ;






 
extern int fputs (const char *__restrict __s, FILE *__restrict __stream);




 
extern int puts (const char *__s);





 
extern int ungetc (int __c, FILE *__stream);





 
extern size_t fread (void *__restrict __ptr, size_t __size,
		     size_t __n, FILE *__restrict __stream) ;



 
extern size_t fwrite (const void *__restrict __ptr, size_t __size,
		      size_t __n, FILE *__restrict __s);








 
extern size_t fread_unlocked (void *__restrict __ptr, size_t __size,
			      size_t __n, FILE *__restrict __stream) ;
extern size_t fwrite_unlocked (const void *__restrict __ptr, size_t __size,
			       size_t __n, FILE *__restrict __stream);






 
extern int fseek (FILE *__stream, long int __off, int __whence);



 
extern long int ftell (FILE *__stream) ;



 
extern void rewind (FILE *__stream);





 




 
extern int fseeko (FILE *__stream, __off_t __off, int __whence);



 
extern __off_t ftello (FILE *__stream) ;





 
extern int fgetpos (FILE *__restrict __stream, fpos_t *__restrict __pos);



 
extern int fsetpos (FILE *__stream, const fpos_t *__pos);




 
extern void clearerr (FILE *__stream) __attribute__ ((__nothrow__ , __leaf__));
 
extern int feof (FILE *__stream) __attribute__ ((__nothrow__ , __leaf__)) ;
 
extern int ferror (FILE *__stream) __attribute__ ((__nothrow__ , __leaf__)) ;


 
extern void clearerr_unlocked (FILE *__stream) __attribute__ ((__nothrow__ , __leaf__));
extern int feof_unlocked (FILE *__stream) __attribute__ ((__nothrow__ , __leaf__)) ;
extern int ferror_unlocked (FILE *__stream) __attribute__ ((__nothrow__ , __leaf__)) ;






 
extern void perror (const char *__s);





 
















 


 

extern int sys_nerr;
extern const char *const sys_errlist[];


 
extern int fileno (FILE *__stream) __attribute__ ((__nothrow__ , __leaf__)) ;

 
extern int fileno_unlocked (FILE *__stream) __attribute__ ((__nothrow__ , __leaf__)) ;





 
extern FILE *popen (const char *__command, const char *__modes) ;




 
extern int pclose (FILE *__stream);


 
extern char *ctermid (char *__s) __attribute__ ((__nothrow__ , __leaf__));






 

 
extern void flockfile (FILE *__stream) __attribute__ ((__nothrow__ , __leaf__));


 
extern int ftrylockfile (FILE *__stream) __attribute__ ((__nothrow__ , __leaf__)) ;

 
extern void funlockfile (FILE *__stream) __attribute__ ((__nothrow__ , __leaf__));



 
















 





 
 
extern __inline __attribute__ ((__gnu_inline__)) int
vprintf (const char *__restrict __fmt, __gnuc_va_list __arg)
{
  return vfprintf (stdout, __fmt, __arg);
}

 
extern __inline __attribute__ ((__gnu_inline__)) int
getchar (void)
{
  return _IO_getc (stdin);
}


 
extern __inline __attribute__ ((__gnu_inline__)) int
fgetc_unlocked (FILE *__fp)
{
  return (__builtin_expect (((__fp)->_IO_read_ptr >= (__fp)->_IO_read_end), 0) ? __uflow (__fp) : *(unsigned char *) (__fp)->_IO_read_ptr++);
}


 
extern __inline __attribute__ ((__gnu_inline__)) int
getc_unlocked (FILE *__fp)
{
  return (__builtin_expect (((__fp)->_IO_read_ptr >= (__fp)->_IO_read_end), 0) ? __uflow (__fp) : *(unsigned char *) (__fp)->_IO_read_ptr++);
}

 
extern __inline __attribute__ ((__gnu_inline__)) int
getchar_unlocked (void)
{
  return (__builtin_expect (((stdin)->_IO_read_ptr >= (stdin)->_IO_read_end), 0) ? __uflow (stdin) : *(unsigned char *) (stdin)->_IO_read_ptr++);
}


 
extern __inline __attribute__ ((__gnu_inline__)) int
putchar (int __c)
{
  return _IO_putc (__c, stdout);
}


 
extern __inline __attribute__ ((__gnu_inline__)) int
fputc_unlocked (int __c, FILE *__stream)
{
  return (__builtin_expect (((__stream)->_IO_write_ptr >= (__stream)->_IO_write_end), 0) ? __overflow (__stream, (unsigned char) (__c)) : (unsigned char) (*(__stream)->_IO_write_ptr++ = (__c)));
}


 
extern __inline __attribute__ ((__gnu_inline__)) int
putc_unlocked (int __c, FILE *__stream)
{
  return (__builtin_expect (((__stream)->_IO_write_ptr >= (__stream)->_IO_write_end), 0) ? __overflow (__stream, (unsigned char) (__c)) : (unsigned char) (*(__stream)->_IO_write_ptr++ = (__c)));
}

 
extern __inline __attribute__ ((__gnu_inline__)) int
putchar_unlocked (int __c)
{
  return (__builtin_expect (((stdout)->_IO_write_ptr >= (stdout)->_IO_write_end), 0) ? __overflow (stdout, (unsigned char) (__c)) : (unsigned char) (*(stdout)->_IO_write_ptr++ = (__c)));
}




 
extern __inline __attribute__ ((__gnu_inline__)) int
__attribute__ ((__nothrow__ , __leaf__)) feof_unlocked (FILE *__stream)
{
  return (((__stream)->_flags & 0x10) != 0);
}

 
extern __inline __attribute__ ((__gnu_inline__)) int
__attribute__ ((__nothrow__ , __leaf__)) ferror_unlocked (FILE *__stream)
{
  return (((__stream)->_flags & 0x20) != 0);
}



 


 



















 



 



 































 



 

 






 




 
















 



 

 

















 



 


 

 

 

 

 

 


 

 

 


















 









 


 
 




 




 
















 



 


 

 

 
















 



 

static __inline unsigned int
__bswap_32 (unsigned int __bsx)
{
  return __builtin_bswap32 (__bsx);
}


 

static __inline __uint64_t
__bswap_64 (__uint64_t __bsx)
{
  return __builtin_bswap64 (__bsx);
}







union wait
  {
    int w_status;
    struct
      {
	unsigned int __w_termsig:7;  
	unsigned int __w_coredump:1;  
	unsigned int __w_retcode:8;  
	unsigned int:16;
      } __wait_terminated;
    struct
      {
	unsigned int __w_stopval:8;  
	unsigned int __w_stopsig:8;  
	unsigned int:16;
      } __wait_stopped;
  };





 





 

 
typedef union
  {
    union wait *__uptr;
    int *__iptr;
  } __WAIT_STATUS __attribute__ ((__transparent_union__));


 


 
typedef struct
  {
    int quot;			 
    int rem;			 
  } div_t;

 
typedef struct
  {
    long int quot;		 
    long int rem;		 
  } ldiv_t;



 
__extension__ typedef struct
  {
    long long int quot;		 
    long long int rem;		 
  } lldiv_t;



 



 


 
extern size_t __ctype_get_mb_cur_max (void) __attribute__ ((__nothrow__ , __leaf__)) ;



 
extern double atof (const char *__nptr)
     __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__pure__)) __attribute__ ((__nonnull__ (1))) ;
 
extern int atoi (const char *__nptr)
     __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__pure__)) __attribute__ ((__nonnull__ (1))) ;
 
extern long int atol (const char *__nptr)
     __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__pure__)) __attribute__ ((__nonnull__ (1))) ;



 
__extension__ extern long long int atoll (const char *__nptr)
     __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__pure__)) __attribute__ ((__nonnull__ (1))) ;



 
extern double strtod (const char *__restrict __nptr,
		      char **__restrict __endptr)
     __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__nonnull__ (1)));



 
extern float strtof (const char *__restrict __nptr,
		     char **__restrict __endptr) __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__nonnull__ (1)));

extern long double strtold (const char *__restrict __nptr,
			    char **__restrict __endptr)
     __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__nonnull__ (1)));



 
extern long int strtol (const char *__restrict __nptr,
			char **__restrict __endptr, int __base)
     __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__nonnull__ (1)));
 
extern unsigned long int strtoul (const char *__restrict __nptr,
				  char **__restrict __endptr, int __base)
     __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__nonnull__ (1)));


 
__extension__
extern long long int strtoq (const char *__restrict __nptr,
			     char **__restrict __endptr, int __base)
     __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__nonnull__ (1)));
 
__extension__
extern unsigned long long int strtouq (const char *__restrict __nptr,
				       char **__restrict __endptr, int __base)
     __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__nonnull__ (1)));


 
__extension__
extern long long int strtoll (const char *__restrict __nptr,
			      char **__restrict __endptr, int __base)
     __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__nonnull__ (1)));
 
__extension__
extern unsigned long long int strtoull (const char *__restrict __nptr,
					char **__restrict __endptr, int __base)
     __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__nonnull__ (1)));






extern __inline __attribute__ ((__gnu_inline__)) int
__attribute__ ((__nothrow__ , __leaf__)) atoi (const char *__nptr)
{
  return (int) strtol (__nptr, (char **) ((void*)0), 10);
}
extern __inline __attribute__ ((__gnu_inline__)) long int
__attribute__ ((__nothrow__ , __leaf__)) atol (const char *__nptr)
{
  return strtol (__nptr, (char **) ((void*)0), 10);
}



__extension__ extern __inline __attribute__ ((__gnu_inline__)) long long int
__attribute__ ((__nothrow__ , __leaf__)) atoll (const char *__nptr)
{
  return strtoll (__nptr, (char **) ((void*)0), 10);
}





 
extern char *l64a (long int __n) __attribute__ ((__nothrow__ , __leaf__)) ;

 
extern long int a64l (const char *__s)
     __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__pure__)) __attribute__ ((__nonnull__ (1))) ;

















 



 






typedef __u_char u_char;
typedef __u_short u_short;
typedef __u_int u_int;
typedef __u_long u_long;
typedef __quad_t quad_t;
typedef __u_quad_t u_quad_t;
typedef __fsid_t fsid_t;

typedef __loff_t loff_t;

typedef __ino_t ino_t;

typedef __dev_t dev_t;

typedef __gid_t gid_t;

typedef __mode_t mode_t;

typedef __nlink_t nlink_t;

typedef __uid_t uid_t;


typedef __pid_t pid_t;

typedef __id_t id_t;


typedef __daddr_t daddr_t;
typedef __caddr_t caddr_t;

typedef __key_t key_t;
















 



 







 
typedef __clock_t clock_t;







 
typedef __time_t time_t;






 
typedef __clockid_t clockid_t;




 
typedef __timer_t timer_t;







































 



 

 






 

 
typedef unsigned long int ulong;
typedef unsigned short int ushort;
typedef unsigned int uint;

 


 

typedef int int8_t __attribute__ ((__mode__ (__QI__)));
typedef int int16_t __attribute__ ((__mode__ (__HI__)));
typedef int int32_t __attribute__ ((__mode__ (__SI__)));
typedef int int64_t __attribute__ ((__mode__ (__DI__)));

typedef unsigned int u_int8_t __attribute__ ((__mode__ (__QI__)));
typedef unsigned int u_int16_t __attribute__ ((__mode__ (__HI__)));
typedef unsigned int u_int32_t __attribute__ ((__mode__ (__SI__)));
typedef unsigned int u_int64_t __attribute__ ((__mode__ (__DI__)));

typedef int register_t __attribute__ ((__mode__ (__word__)));



 


 

 
















 

 



 

 















 


 


 







 

















 


typedef int __sig_atomic_t;

 

typedef struct
  {
    unsigned long int __val[(1024 / (8 * sizeof (unsigned long int)))];
  } __sigset_t;







 


typedef __sigset_t sigset_t;

 















 



 












 
struct timespec
  {
    __time_t tv_sec;		 
    __syscall_slong_t tv_nsec;	 
  };




















 



 



 
struct timeval
  {
    __time_t tv_sec;		 
    __suseconds_t tv_usec;	 
  };



typedef __suseconds_t suseconds_t;


 
typedef long int __fd_mask;

 
 

 
typedef struct
  {
    
 
    __fd_mask __fds_bits[1024 / (8 * (int) sizeof (__fd_mask))];
  } fd_set;

 

 
typedef __fd_mask fd_mask;

 


 











 
extern int select (int __nfds, fd_set *__restrict __readfds,
		   fd_set *__restrict __writefds,
		   fd_set *__restrict __exceptfds,
		   struct timeval *__restrict __timeout);






 
extern int pselect (int __nfds, fd_set *__restrict __readfds,
		    fd_set *__restrict __writefds,
		    fd_set *__restrict __exceptfds,
		    const struct timespec *__restrict __timeout,
		    const __sigset_t *__restrict __sigmask);


 




 

















 





 


__extension__
extern unsigned int gnu_dev_major (unsigned long long int __dev)
     __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__const__));
__extension__
extern unsigned int gnu_dev_minor (unsigned long long int __dev)
     __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__const__));
__extension__
extern unsigned long long int gnu_dev_makedev (unsigned int __major,
					       unsigned int __minor)
     __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__const__));

__extension__ extern __inline __attribute__ ((__gnu_inline__)) __attribute__ ((__const__)) unsigned int
__attribute__ ((__nothrow__ , __leaf__)) gnu_dev_major (unsigned long long int __dev)
{
  return ((__dev >> 8) & 0xfff) | ((unsigned int) (__dev >> 32) & ~0xfff);
}

__extension__ extern __inline __attribute__ ((__gnu_inline__)) __attribute__ ((__const__)) unsigned int
__attribute__ ((__nothrow__ , __leaf__)) gnu_dev_minor (unsigned long long int __dev)
{
  return (__dev & 0xff) | ((unsigned int) (__dev >> 12) & ~0xff);
}

__extension__ extern __inline __attribute__ ((__gnu_inline__)) __attribute__ ((__const__)) unsigned long long int
__attribute__ ((__nothrow__ , __leaf__)) gnu_dev_makedev (unsigned int __major, unsigned int __minor)
{
  return ((__minor & 0xff) | ((__major & 0xfff) << 8)
	  | (((unsigned long long int) (__minor & ~0xff)) << 12)
	  | (((unsigned long long int) (__major & ~0xfff)) << 32));
}


 



typedef __blksize_t blksize_t;

 
typedef __blkcnt_t blkcnt_t;	  
typedef __fsblkcnt_t fsblkcnt_t;  
typedef __fsfilcnt_t fsfilcnt_t;  



 















 


 


 




 
typedef unsigned long int pthread_t;


union pthread_attr_t
{
  char __size[56];
  long int __align;
};
typedef union pthread_attr_t pthread_attr_t;


typedef struct __pthread_internal_list
{
  struct __pthread_internal_list *__prev;
  struct __pthread_internal_list *__next;
} __pthread_list_t;



 
typedef union
{
  struct __pthread_mutex_s
  {
    int __lock;
    unsigned int __count;
    int __owner;
    unsigned int __nusers;
    
 
    int __kind;
    short __spins;
    short __elision;
    __pthread_list_t __list;
 
  } __data;
  char __size[40];
  long int __align;
} pthread_mutex_t;

typedef union
{
  char __size[4];
  int __align;
} pthread_mutexattr_t;



 
typedef union
{
  struct
  {
    int __lock;
    unsigned int __futex;
    __extension__ unsigned long long int __total_seq;
    __extension__ unsigned long long int __wakeup_seq;
    __extension__ unsigned long long int __woken_seq;
    void *__mutex;
    unsigned int __nwaiters;
    unsigned int __broadcast_seq;
  } __data;
  char __size[48];
  __extension__ long long int __align;
} pthread_cond_t;

typedef union
{
  char __size[4];
  int __align;
} pthread_condattr_t;


 
typedef unsigned int pthread_key_t;


 
typedef int pthread_once_t;



 
typedef union
{
  struct
  {
    int __lock;
    unsigned int __nr_readers;
    unsigned int __readers_wakeup;
    unsigned int __writer_wakeup;
    unsigned int __nr_readers_queued;
    unsigned int __nr_writers_queued;
    int __writer;
    int __shared;
    unsigned long int __pad1;
    unsigned long int __pad2;
    
 
    unsigned int __flags;
  } __data;
  char __size[56];
  long int __align;
} pthread_rwlock_t;

typedef union
{
  char __size[8];
  long int __align;
} pthread_rwlockattr_t;


 
typedef volatile int pthread_spinlock_t;



 
typedef union
{
  char __size[32];
  long int __align;
} pthread_barrier_t;

typedef union
{
  char __size[4];
  int __align;
} pthread_barrierattr_t;










 
 
extern long int random (void) __attribute__ ((__nothrow__ , __leaf__));

 
extern void srandom (unsigned int __seed) __attribute__ ((__nothrow__ , __leaf__));




 
extern char *initstate (unsigned int __seed, char *__statebuf,
			size_t __statelen) __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__nonnull__ (2)));


 
extern char *setstate (char *__statebuf) __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__nonnull__ (1)));




 

struct random_data
  {
    int32_t *fptr;		 
    int32_t *rptr;		 
    int32_t *state;		 
    int rand_type;		 
    int rand_deg;		 
    int rand_sep;		 
    int32_t *end_ptr;		 
  };

extern int random_r (struct random_data *__restrict __buf,
		     int32_t *__restrict __result) __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__nonnull__ (1, 2)));

extern int srandom_r (unsigned int __seed, struct random_data *__buf)
     __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__nonnull__ (2)));

extern int initstate_r (unsigned int __seed, char *__restrict __statebuf,
			size_t __statelen,
			struct random_data *__restrict __buf)
     __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__nonnull__ (2, 4)));

extern int setstate_r (char *__restrict __statebuf,
		       struct random_data *__restrict __buf)
     __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__nonnull__ (1, 2)));



 
extern int rand (void) __attribute__ ((__nothrow__ , __leaf__));
 
extern void srand (unsigned int __seed) __attribute__ ((__nothrow__ , __leaf__));


 
extern int rand_r (unsigned int *__seed) __attribute__ ((__nothrow__ , __leaf__));


 

 
extern double drand48 (void) __attribute__ ((__nothrow__ , __leaf__));
extern double erand48 (unsigned short int __xsubi[3]) __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__nonnull__ (1)));

 
extern long int lrand48 (void) __attribute__ ((__nothrow__ , __leaf__));
extern long int nrand48 (unsigned short int __xsubi[3])
     __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__nonnull__ (1)));

 
extern long int mrand48 (void) __attribute__ ((__nothrow__ , __leaf__));
extern long int jrand48 (unsigned short int __xsubi[3])
     __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__nonnull__ (1)));

 
extern void srand48 (long int __seedval) __attribute__ ((__nothrow__ , __leaf__));
extern unsigned short int *seed48 (unsigned short int __seed16v[3])
     __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__nonnull__ (1)));
extern void lcong48 (unsigned short int __param[7]) __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__nonnull__ (1)));



 
struct drand48_data
  {
    unsigned short int __x[3];	 
    unsigned short int __old_x[3];  
    unsigned short int __c;	 
    unsigned short int __init;	 
    unsigned long long int __a;	 
  };

 
extern int drand48_r (struct drand48_data *__restrict __buffer,
		      double *__restrict __result) __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__nonnull__ (1, 2)));
extern int erand48_r (unsigned short int __xsubi[3],
		      struct drand48_data *__restrict __buffer,
		      double *__restrict __result) __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__nonnull__ (1, 2)));

 
extern int lrand48_r (struct drand48_data *__restrict __buffer,
		      long int *__restrict __result)
     __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__nonnull__ (1, 2)));
extern int nrand48_r (unsigned short int __xsubi[3],
		      struct drand48_data *__restrict __buffer,
		      long int *__restrict __result)
     __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__nonnull__ (1, 2)));

 
extern int mrand48_r (struct drand48_data *__restrict __buffer,
		      long int *__restrict __result)
     __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__nonnull__ (1, 2)));
extern int jrand48_r (unsigned short int __xsubi[3],
		      struct drand48_data *__restrict __buffer,
		      long int *__restrict __result)
     __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__nonnull__ (1, 2)));

 
extern int srand48_r (long int __seedval, struct drand48_data *__buffer)
     __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__nonnull__ (2)));

extern int seed48_r (unsigned short int __seed16v[3],
		     struct drand48_data *__buffer) __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__nonnull__ (1, 2)));

extern int lcong48_r (unsigned short int __param[7],
		      struct drand48_data *__buffer)
     __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__nonnull__ (1, 2)));



 
extern void *malloc (size_t __size) __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__malloc__)) ;
 
extern void *calloc (size_t __nmemb, size_t __size)
     __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__malloc__)) ;




 


 
extern void *realloc (void *__ptr, size_t __size)
     __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__warn_unused_result__));
 
extern void free (void *__ptr) __attribute__ ((__nothrow__ , __leaf__));


 
extern void cfree (void *__ptr) __attribute__ ((__nothrow__ , __leaf__));
















 


































 



 

 






 



 

 
extern void *alloca (size_t __size) __attribute__ ((__nothrow__ , __leaf__));





 
extern void *valloc (size_t __size) __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__malloc__)) ;

 
extern int posix_memalign (void **__memptr, size_t __alignment, size_t __size)
     __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__nonnull__ (1))) ;

 
extern void *aligned_alloc (size_t __alignment, size_t __size)
     __attribute__ ((__nothrow__ , __leaf__))  __attribute__ ((__malloc__, __alloc_size__ (2)));


 
extern void abort (void) __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__noreturn__));


 
extern int atexit (void (*__func) (void)) __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__nonnull__ (1)));

 
extern int at_quick_exit (void (*__func) (void)) __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__nonnull__ (1)));



 
extern int on_exit (void (*__func) (int __status, void *__arg), void *__arg)
     __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__nonnull__ (1)));




 
extern void exit (int __status) __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__noreturn__));



 
extern void quick_exit (int __status) __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__noreturn__));




 
extern void _Exit (int __status) __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__noreturn__));




 
extern char *getenv (const char *__name) __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__nonnull__ (1))) ;



 

 
extern int putenv (char *__string) __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__nonnull__ (1)));


 
extern int setenv (const char *__name, const char *__value, int __replace)
     __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__nonnull__ (2)));

 
extern int unsetenv (const char *__name) __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__nonnull__ (1)));



 
extern int clearenv (void) __attribute__ ((__nothrow__ , __leaf__));






 
extern char *mktemp (char *__template) __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__nonnull__ (1)));








 
extern int mkstemp (char *__template) __attribute__ ((__nonnull__ (1))) ;






 
extern int mkstemps (char *__template, int __suffixlen) __attribute__ ((__nonnull__ (1))) ;





 
extern char *mkdtemp (char *__template) __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__nonnull__ (1))) ;







 
extern int system (const char *__command) ;








 
extern char *realpath (const char *__restrict __name,
		       char *__restrict __resolved) __attribute__ ((__nothrow__ , __leaf__)) ;


 
typedef int (*__compar_fn_t) (const void *, const void *);




 
extern void *bsearch (const void *__key, const void *__base,
		      size_t __nmemb, size_t __size, __compar_fn_t __compar)
     __attribute__ ((__nonnull__ (1, 2, 5))) ;


 
extern void qsort (void *__base, size_t __nmemb, size_t __size,
		   __compar_fn_t __compar) __attribute__ ((__nonnull__ (1, 4)));


 
extern int abs (int __x) __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__const__)) ;
extern long int labs (long int __x) __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__const__)) ;


__extension__ extern long long int llabs (long long int __x)
     __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__const__)) ;




 
 
extern div_t div (int __numer, int __denom)
     __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__const__)) ;
extern ldiv_t ldiv (long int __numer, long int __denom)
     __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__const__)) ;



__extension__ extern lldiv_t lldiv (long long int __numer,
				    long long int __denom)
     __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__const__)) ;




 



 
extern char *ecvt (double __value, int __ndigit, int *__restrict __decpt,
		   int *__restrict __sign) __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__nonnull__ (3, 4))) ;



 
extern char *fcvt (double __value, int __ndigit, int *__restrict __decpt,
		   int *__restrict __sign) __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__nonnull__ (3, 4))) ;



 
extern char *gcvt (double __value, int __ndigit, char *__buf)
     __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__nonnull__ (3))) ;


 
extern char *qecvt (long double __value, int __ndigit,
		    int *__restrict __decpt, int *__restrict __sign)
     __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__nonnull__ (3, 4))) ;
extern char *qfcvt (long double __value, int __ndigit,
		    int *__restrict __decpt, int *__restrict __sign)
     __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__nonnull__ (3, 4))) ;
extern char *qgcvt (long double __value, int __ndigit, char *__buf)
     __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__nonnull__ (3))) ;



 
extern int ecvt_r (double __value, int __ndigit, int *__restrict __decpt,
		   int *__restrict __sign, char *__restrict __buf,
		   size_t __len) __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__nonnull__ (3, 4, 5)));
extern int fcvt_r (double __value, int __ndigit, int *__restrict __decpt,
		   int *__restrict __sign, char *__restrict __buf,
		   size_t __len) __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__nonnull__ (3, 4, 5)));

extern int qecvt_r (long double __value, int __ndigit,
		    int *__restrict __decpt, int *__restrict __sign,
		    char *__restrict __buf, size_t __len)
     __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__nonnull__ (3, 4, 5)));
extern int qfcvt_r (long double __value, int __ndigit,
		    int *__restrict __decpt, int *__restrict __sign,
		    char *__restrict __buf, size_t __len)
     __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__nonnull__ (3, 4, 5)));




 
extern int mblen (const char *__s, size_t __n) __attribute__ ((__nothrow__ , __leaf__)) ;

 
extern int mbtowc (wchar_t *__restrict __pwc,
		   const char *__restrict __s, size_t __n) __attribute__ ((__nothrow__ , __leaf__)) ;

 
extern int wctomb (char *__s, wchar_t __wchar) __attribute__ ((__nothrow__ , __leaf__)) ;


 
extern size_t mbstowcs (wchar_t *__restrict  __pwcs,
			const char *__restrict __s, size_t __n) __attribute__ ((__nothrow__ , __leaf__));
 
extern size_t wcstombs (char *__restrict __s,
			const wchar_t *__restrict __pwcs, size_t __n)
     __attribute__ ((__nothrow__ , __leaf__));






 
extern int rpmatch (const char *__response) __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__nonnull__ (1))) ;







 
extern int getsubopt (char **__restrict __optionp,
		      char *const *__restrict __tokens,
		      char **__restrict __valuep)
     __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__nonnull__ (1, 2, 3))) ;




 






 
extern int getloadavg (double __loadavg[], int __nelem)
     __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__nonnull__ (1)));

















 



extern __inline __attribute__ ((__gnu_inline__)) double
__attribute__ ((__nothrow__ , __leaf__)) atof (const char *__nptr)
{
  return strtod (__nptr, (char **) ((void*)0));
}


 




 








 







 

 

















 



 






 


















 


 



















 


 



















 



 
















 


 


 
















 



 




 















 




 
typedef float float_t;		 
typedef double double_t;	
 

 


 







 


















 

























 



 


 
extern double acos (double __x) __attribute__ ((__nothrow__ , __leaf__)); extern double __acos (double __x) __attribute__ ((__nothrow__ , __leaf__));
 
extern double asin (double __x) __attribute__ ((__nothrow__ , __leaf__)); extern double __asin (double __x) __attribute__ ((__nothrow__ , __leaf__));
 
extern double atan (double __x) __attribute__ ((__nothrow__ , __leaf__)); extern double __atan (double __x) __attribute__ ((__nothrow__ , __leaf__));
 
extern double atan2 (double __y, double __x) __attribute__ ((__nothrow__ , __leaf__)); extern double __atan2 (double __y, double __x) __attribute__ ((__nothrow__ , __leaf__));

 
extern double cos (double __x) __attribute__ ((__nothrow__ , __leaf__)); extern double __cos (double __x) __attribute__ ((__nothrow__ , __leaf__));
 
extern double sin (double __x) __attribute__ ((__nothrow__ , __leaf__)); extern double __sin (double __x) __attribute__ ((__nothrow__ , __leaf__));
 
extern double tan (double __x) __attribute__ ((__nothrow__ , __leaf__)); extern double __tan (double __x) __attribute__ ((__nothrow__ , __leaf__));

 

 
extern double cosh (double __x) __attribute__ ((__nothrow__ , __leaf__)); extern double __cosh (double __x) __attribute__ ((__nothrow__ , __leaf__));
 
extern double sinh (double __x) __attribute__ ((__nothrow__ , __leaf__)); extern double __sinh (double __x) __attribute__ ((__nothrow__ , __leaf__));
 
extern double tanh (double __x) __attribute__ ((__nothrow__ , __leaf__)); extern double __tanh (double __x) __attribute__ ((__nothrow__ , __leaf__));




 
extern double acosh (double __x) __attribute__ ((__nothrow__ , __leaf__)); extern double __acosh (double __x) __attribute__ ((__nothrow__ , __leaf__));
 
extern double asinh (double __x) __attribute__ ((__nothrow__ , __leaf__)); extern double __asinh (double __x) __attribute__ ((__nothrow__ , __leaf__));
 
extern double atanh (double __x) __attribute__ ((__nothrow__ , __leaf__)); extern double __atanh (double __x) __attribute__ ((__nothrow__ , __leaf__));


 


 
extern double exp (double __x) __attribute__ ((__nothrow__ , __leaf__)); extern double __exp (double __x) __attribute__ ((__nothrow__ , __leaf__));

 
extern double frexp (double __x, int *__exponent) __attribute__ ((__nothrow__ , __leaf__)); extern double __frexp (double __x, int *__exponent) __attribute__ ((__nothrow__ , __leaf__));

 
extern double ldexp (double __x, int __exponent) __attribute__ ((__nothrow__ , __leaf__)); extern double __ldexp (double __x, int __exponent) __attribute__ ((__nothrow__ , __leaf__));

 
extern double log (double __x) __attribute__ ((__nothrow__ , __leaf__)); extern double __log (double __x) __attribute__ ((__nothrow__ , __leaf__));

 
extern double log10 (double __x) __attribute__ ((__nothrow__ , __leaf__)); extern double __log10 (double __x) __attribute__ ((__nothrow__ , __leaf__));

 
extern double modf (double __x, double *__iptr) __attribute__ ((__nothrow__ , __leaf__)); extern double __modf (double __x, double *__iptr) __attribute__ ((__nothrow__ , __leaf__))
     __attribute__ ((__nonnull__ (2)));




 
extern double expm1 (double __x) __attribute__ ((__nothrow__ , __leaf__)); extern double __expm1 (double __x) __attribute__ ((__nothrow__ , __leaf__));

 
extern double log1p (double __x) __attribute__ ((__nothrow__ , __leaf__)); extern double __log1p (double __x) __attribute__ ((__nothrow__ , __leaf__));

 
extern double logb (double __x) __attribute__ ((__nothrow__ , __leaf__)); extern double __logb (double __x) __attribute__ ((__nothrow__ , __leaf__));



 
extern double exp2 (double __x) __attribute__ ((__nothrow__ , __leaf__)); extern double __exp2 (double __x) __attribute__ ((__nothrow__ , __leaf__));

 
extern double log2 (double __x) __attribute__ ((__nothrow__ , __leaf__)); extern double __log2 (double __x) __attribute__ ((__nothrow__ , __leaf__));



 


 
extern double pow (double __x, double __y) __attribute__ ((__nothrow__ , __leaf__)); extern double __pow (double __x, double __y) __attribute__ ((__nothrow__ , __leaf__));

 
extern double sqrt (double __x) __attribute__ ((__nothrow__ , __leaf__)); extern double __sqrt (double __x) __attribute__ ((__nothrow__ , __leaf__));



 
extern double hypot (double __x, double __y) __attribute__ ((__nothrow__ , __leaf__)); extern double __hypot (double __x, double __y) __attribute__ ((__nothrow__ , __leaf__));



 
extern double cbrt (double __x) __attribute__ ((__nothrow__ , __leaf__)); extern double __cbrt (double __x) __attribute__ ((__nothrow__ , __leaf__));



 


 
extern double ceil (double __x) __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__const__)); extern double __ceil (double __x) __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__const__));

 
extern double fabs (double __x) __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__const__)); extern double __fabs (double __x) __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__const__));

 
extern double floor (double __x) __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__const__)); extern double __floor (double __x) __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__const__));

 
extern double fmod (double __x, double __y) __attribute__ ((__nothrow__ , __leaf__)); extern double __fmod (double __x, double __y) __attribute__ ((__nothrow__ , __leaf__));



 
extern int __isinf (double __value) __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__const__));

 
extern int __finite (double __value) __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__const__));



 
extern int isinf (double __value) __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__const__));

 
extern int finite (double __value) __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__const__));

 
extern double drem (double __x, double __y) __attribute__ ((__nothrow__ , __leaf__)); extern double __drem (double __x, double __y) __attribute__ ((__nothrow__ , __leaf__));


 
extern double significand (double __x) __attribute__ ((__nothrow__ , __leaf__)); extern double __significand (double __x) __attribute__ ((__nothrow__ , __leaf__));


 
extern double copysign (double __x, double __y) __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__const__)); extern double __copysign (double __x, double __y) __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__const__));



 
extern double nan (const char *__tagb) __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__const__)); extern double __nan (const char *__tagb) __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__const__));



 
extern int __isnan (double __value) __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__const__));

 
extern int isnan (double __value) __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__const__));

 
extern double j0 (double) __attribute__ ((__nothrow__ , __leaf__)); extern double __j0 (double) __attribute__ ((__nothrow__ , __leaf__));
extern double j1 (double) __attribute__ ((__nothrow__ , __leaf__)); extern double __j1 (double) __attribute__ ((__nothrow__ , __leaf__));
extern double jn (int, double) __attribute__ ((__nothrow__ , __leaf__)); extern double __jn (int, double) __attribute__ ((__nothrow__ , __leaf__));
extern double y0 (double) __attribute__ ((__nothrow__ , __leaf__)); extern double __y0 (double) __attribute__ ((__nothrow__ , __leaf__));
extern double y1 (double) __attribute__ ((__nothrow__ , __leaf__)); extern double __y1 (double) __attribute__ ((__nothrow__ , __leaf__));
extern double yn (int, double) __attribute__ ((__nothrow__ , __leaf__)); extern double __yn (int, double) __attribute__ ((__nothrow__ , __leaf__));



 
extern double erf (double) __attribute__ ((__nothrow__ , __leaf__)); extern double __erf (double) __attribute__ ((__nothrow__ , __leaf__));
extern double erfc (double) __attribute__ ((__nothrow__ , __leaf__)); extern double __erfc (double) __attribute__ ((__nothrow__ , __leaf__));
extern double lgamma (double) __attribute__ ((__nothrow__ , __leaf__)); extern double __lgamma (double) __attribute__ ((__nothrow__ , __leaf__));



 
extern double tgamma (double) __attribute__ ((__nothrow__ , __leaf__)); extern double __tgamma (double) __attribute__ ((__nothrow__ , __leaf__));


 
extern double gamma (double) __attribute__ ((__nothrow__ , __leaf__)); extern double __gamma (double) __attribute__ ((__nothrow__ , __leaf__));



 
extern double lgamma_r (double, int *__signgamp) __attribute__ ((__nothrow__ , __leaf__)); extern double __lgamma_r (double, int *__signgamp) __attribute__ ((__nothrow__ , __leaf__));




 
extern double rint (double __x) __attribute__ ((__nothrow__ , __leaf__)); extern double __rint (double __x) __attribute__ ((__nothrow__ , __leaf__));

 
extern double nextafter (double __x, double __y) __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__const__)); extern double __nextafter (double __x, double __y) __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__const__));
extern double nexttoward (double __x, long double __y) __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__const__)); extern double __nexttoward (double __x, long double __y) __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__const__));

 
extern double remainder (double __x, double __y) __attribute__ ((__nothrow__ , __leaf__)); extern double __remainder (double __x, double __y) __attribute__ ((__nothrow__ , __leaf__));

 
extern double scalbn (double __x, int __n) __attribute__ ((__nothrow__ , __leaf__)); extern double __scalbn (double __x, int __n) __attribute__ ((__nothrow__ , __leaf__));

 
extern int ilogb (double __x) __attribute__ ((__nothrow__ , __leaf__)); extern int __ilogb (double __x) __attribute__ ((__nothrow__ , __leaf__));

 
extern double scalbln (double __x, long int __n) __attribute__ ((__nothrow__ , __leaf__)); extern double __scalbln (double __x, long int __n) __attribute__ ((__nothrow__ , __leaf__));


 
extern double nearbyint (double __x) __attribute__ ((__nothrow__ , __leaf__)); extern double __nearbyint (double __x) __attribute__ ((__nothrow__ , __leaf__));


 
extern double round (double __x) __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__const__)); extern double __round (double __x) __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__const__));


 
extern double trunc (double __x) __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__const__)); extern double __trunc (double __x) __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__const__));



 
extern double remquo (double __x, double __y, int *__quo) __attribute__ ((__nothrow__ , __leaf__)); extern double __remquo (double __x, double __y, int *__quo) __attribute__ ((__nothrow__ , __leaf__));


 


 
extern long int lrint (double __x) __attribute__ ((__nothrow__ , __leaf__)); extern long int __lrint (double __x) __attribute__ ((__nothrow__ , __leaf__));
extern long long int llrint (double __x) __attribute__ ((__nothrow__ , __leaf__)); extern long long int __llrint (double __x) __attribute__ ((__nothrow__ , __leaf__));


 
extern long int lround (double __x) __attribute__ ((__nothrow__ , __leaf__)); extern long int __lround (double __x) __attribute__ ((__nothrow__ , __leaf__));
extern long long int llround (double __x) __attribute__ ((__nothrow__ , __leaf__)); extern long long int __llround (double __x) __attribute__ ((__nothrow__ , __leaf__));


 
extern double fdim (double __x, double __y) __attribute__ ((__nothrow__ , __leaf__)); extern double __fdim (double __x, double __y) __attribute__ ((__nothrow__ , __leaf__));

 
extern double fmax (double __x, double __y) __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__const__)); extern double __fmax (double __x, double __y) __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__const__));

 
extern double fmin (double __x, double __y) __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__const__)); extern double __fmin (double __x, double __y) __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__const__));


 
extern int __fpclassify (double __value) __attribute__ ((__nothrow__ , __leaf__))
     __attribute__ ((__const__));

 
extern int __signbit (double __value) __attribute__ ((__nothrow__ , __leaf__))
     __attribute__ ((__const__));


 
extern double fma (double __x, double __y, double __z) __attribute__ ((__nothrow__ , __leaf__)); extern double __fma (double __x, double __y, double __z) __attribute__ ((__nothrow__ , __leaf__));



 
extern double scalb (double __x, double __n) __attribute__ ((__nothrow__ , __leaf__)); extern double __scalb (double __x, double __n) __attribute__ ((__nothrow__ , __leaf__));




 

















 

























 



 


 
extern float acosf (float __x) __attribute__ ((__nothrow__ , __leaf__)); extern float __acosf (float __x) __attribute__ ((__nothrow__ , __leaf__));
 
extern float asinf (float __x) __attribute__ ((__nothrow__ , __leaf__)); extern float __asinf (float __x) __attribute__ ((__nothrow__ , __leaf__));
 
extern float atanf (float __x) __attribute__ ((__nothrow__ , __leaf__)); extern float __atanf (float __x) __attribute__ ((__nothrow__ , __leaf__));
 
extern float atan2f (float __y, float __x) __attribute__ ((__nothrow__ , __leaf__)); extern float __atan2f (float __y, float __x) __attribute__ ((__nothrow__ , __leaf__));

 
extern float cosf (float __x) __attribute__ ((__nothrow__ , __leaf__)); extern float __cosf (float __x) __attribute__ ((__nothrow__ , __leaf__));
 
extern float sinf (float __x) __attribute__ ((__nothrow__ , __leaf__)); extern float __sinf (float __x) __attribute__ ((__nothrow__ , __leaf__));
 
extern float tanf (float __x) __attribute__ ((__nothrow__ , __leaf__)); extern float __tanf (float __x) __attribute__ ((__nothrow__ , __leaf__));

 

 
extern float coshf (float __x) __attribute__ ((__nothrow__ , __leaf__)); extern float __coshf (float __x) __attribute__ ((__nothrow__ , __leaf__));
 
extern float sinhf (float __x) __attribute__ ((__nothrow__ , __leaf__)); extern float __sinhf (float __x) __attribute__ ((__nothrow__ , __leaf__));
 
extern float tanhf (float __x) __attribute__ ((__nothrow__ , __leaf__)); extern float __tanhf (float __x) __attribute__ ((__nothrow__ , __leaf__));




 
extern float acoshf (float __x) __attribute__ ((__nothrow__ , __leaf__)); extern float __acoshf (float __x) __attribute__ ((__nothrow__ , __leaf__));
 
extern float asinhf (float __x) __attribute__ ((__nothrow__ , __leaf__)); extern float __asinhf (float __x) __attribute__ ((__nothrow__ , __leaf__));
 
extern float atanhf (float __x) __attribute__ ((__nothrow__ , __leaf__)); extern float __atanhf (float __x) __attribute__ ((__nothrow__ , __leaf__));


 


 
extern float expf (float __x) __attribute__ ((__nothrow__ , __leaf__)); extern float __expf (float __x) __attribute__ ((__nothrow__ , __leaf__));

 
extern float frexpf (float __x, int *__exponent) __attribute__ ((__nothrow__ , __leaf__)); extern float __frexpf (float __x, int *__exponent) __attribute__ ((__nothrow__ , __leaf__));

 
extern float ldexpf (float __x, int __exponent) __attribute__ ((__nothrow__ , __leaf__)); extern float __ldexpf (float __x, int __exponent) __attribute__ ((__nothrow__ , __leaf__));

 
extern float logf (float __x) __attribute__ ((__nothrow__ , __leaf__)); extern float __logf (float __x) __attribute__ ((__nothrow__ , __leaf__));

 
extern float log10f (float __x) __attribute__ ((__nothrow__ , __leaf__)); extern float __log10f (float __x) __attribute__ ((__nothrow__ , __leaf__));

 
extern float modff (float __x, float *__iptr) __attribute__ ((__nothrow__ , __leaf__)); extern float __modff (float __x, float *__iptr) __attribute__ ((__nothrow__ , __leaf__))
     __attribute__ ((__nonnull__ (2)));




 
extern float expm1f (float __x) __attribute__ ((__nothrow__ , __leaf__)); extern float __expm1f (float __x) __attribute__ ((__nothrow__ , __leaf__));

 
extern float log1pf (float __x) __attribute__ ((__nothrow__ , __leaf__)); extern float __log1pf (float __x) __attribute__ ((__nothrow__ , __leaf__));

 
extern float logbf (float __x) __attribute__ ((__nothrow__ , __leaf__)); extern float __logbf (float __x) __attribute__ ((__nothrow__ , __leaf__));



 
extern float exp2f (float __x) __attribute__ ((__nothrow__ , __leaf__)); extern float __exp2f (float __x) __attribute__ ((__nothrow__ , __leaf__));

 
extern float log2f (float __x) __attribute__ ((__nothrow__ , __leaf__)); extern float __log2f (float __x) __attribute__ ((__nothrow__ , __leaf__));



 


 
extern float powf (float __x, float __y) __attribute__ ((__nothrow__ , __leaf__)); extern float __powf (float __x, float __y) __attribute__ ((__nothrow__ , __leaf__));

 
extern float sqrtf (float __x) __attribute__ ((__nothrow__ , __leaf__)); extern float __sqrtf (float __x) __attribute__ ((__nothrow__ , __leaf__));



 
extern float hypotf (float __x, float __y) __attribute__ ((__nothrow__ , __leaf__)); extern float __hypotf (float __x, float __y) __attribute__ ((__nothrow__ , __leaf__));



 
extern float cbrtf (float __x) __attribute__ ((__nothrow__ , __leaf__)); extern float __cbrtf (float __x) __attribute__ ((__nothrow__ , __leaf__));



 


 
extern float ceilf (float __x) __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__const__)); extern float __ceilf (float __x) __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__const__));

 
extern float fabsf (float __x) __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__const__)); extern float __fabsf (float __x) __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__const__));

 
extern float floorf (float __x) __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__const__)); extern float __floorf (float __x) __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__const__));

 
extern float fmodf (float __x, float __y) __attribute__ ((__nothrow__ , __leaf__)); extern float __fmodf (float __x, float __y) __attribute__ ((__nothrow__ , __leaf__));



 
extern int __isinff (float __value) __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__const__));

 
extern int __finitef (float __value) __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__const__));



 
extern int isinff (float __value) __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__const__));

 
extern int finitef (float __value) __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__const__));

 
extern float dremf (float __x, float __y) __attribute__ ((__nothrow__ , __leaf__)); extern float __dremf (float __x, float __y) __attribute__ ((__nothrow__ , __leaf__));


 
extern float significandf (float __x) __attribute__ ((__nothrow__ , __leaf__)); extern float __significandf (float __x) __attribute__ ((__nothrow__ , __leaf__));


 
extern float copysignf (float __x, float __y) __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__const__)); extern float __copysignf (float __x, float __y) __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__const__));



 
extern float nanf (const char *__tagb) __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__const__)); extern float __nanf (const char *__tagb) __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__const__));



 
extern int __isnanf (float __value) __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__const__));

 
extern int isnanf (float __value) __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__const__));

 
extern float j0f (float) __attribute__ ((__nothrow__ , __leaf__)); extern float __j0f (float) __attribute__ ((__nothrow__ , __leaf__));
extern float j1f (float) __attribute__ ((__nothrow__ , __leaf__)); extern float __j1f (float) __attribute__ ((__nothrow__ , __leaf__));
extern float jnf (int, float) __attribute__ ((__nothrow__ , __leaf__)); extern float __jnf (int, float) __attribute__ ((__nothrow__ , __leaf__));
extern float y0f (float) __attribute__ ((__nothrow__ , __leaf__)); extern float __y0f (float) __attribute__ ((__nothrow__ , __leaf__));
extern float y1f (float) __attribute__ ((__nothrow__ , __leaf__)); extern float __y1f (float) __attribute__ ((__nothrow__ , __leaf__));
extern float ynf (int, float) __attribute__ ((__nothrow__ , __leaf__)); extern float __ynf (int, float) __attribute__ ((__nothrow__ , __leaf__));



 
extern float erff (float) __attribute__ ((__nothrow__ , __leaf__)); extern float __erff (float) __attribute__ ((__nothrow__ , __leaf__));
extern float erfcf (float) __attribute__ ((__nothrow__ , __leaf__)); extern float __erfcf (float) __attribute__ ((__nothrow__ , __leaf__));
extern float lgammaf (float) __attribute__ ((__nothrow__ , __leaf__)); extern float __lgammaf (float) __attribute__ ((__nothrow__ , __leaf__));



 
extern float tgammaf (float) __attribute__ ((__nothrow__ , __leaf__)); extern float __tgammaf (float) __attribute__ ((__nothrow__ , __leaf__));


 
extern float gammaf (float) __attribute__ ((__nothrow__ , __leaf__)); extern float __gammaf (float) __attribute__ ((__nothrow__ , __leaf__));



 
extern float lgammaf_r (float, int *__signgamp) __attribute__ ((__nothrow__ , __leaf__)); extern float __lgammaf_r (float, int *__signgamp) __attribute__ ((__nothrow__ , __leaf__));




 
extern float rintf (float __x) __attribute__ ((__nothrow__ , __leaf__)); extern float __rintf (float __x) __attribute__ ((__nothrow__ , __leaf__));

 
extern float nextafterf (float __x, float __y) __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__const__)); extern float __nextafterf (float __x, float __y) __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__const__));
extern float nexttowardf (float __x, long double __y) __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__const__)); extern float __nexttowardf (float __x, long double __y) __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__const__));

 
extern float remainderf (float __x, float __y) __attribute__ ((__nothrow__ , __leaf__)); extern float __remainderf (float __x, float __y) __attribute__ ((__nothrow__ , __leaf__));

 
extern float scalbnf (float __x, int __n) __attribute__ ((__nothrow__ , __leaf__)); extern float __scalbnf (float __x, int __n) __attribute__ ((__nothrow__ , __leaf__));

 
extern int ilogbf (float __x) __attribute__ ((__nothrow__ , __leaf__)); extern int __ilogbf (float __x) __attribute__ ((__nothrow__ , __leaf__));

 
extern float scalblnf (float __x, long int __n) __attribute__ ((__nothrow__ , __leaf__)); extern float __scalblnf (float __x, long int __n) __attribute__ ((__nothrow__ , __leaf__));


 
extern float nearbyintf (float __x) __attribute__ ((__nothrow__ , __leaf__)); extern float __nearbyintf (float __x) __attribute__ ((__nothrow__ , __leaf__));


 
extern float roundf (float __x) __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__const__)); extern float __roundf (float __x) __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__const__));


 
extern float truncf (float __x) __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__const__)); extern float __truncf (float __x) __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__const__));



 
extern float remquof (float __x, float __y, int *__quo) __attribute__ ((__nothrow__ , __leaf__)); extern float __remquof (float __x, float __y, int *__quo) __attribute__ ((__nothrow__ , __leaf__));


 


 
extern long int lrintf (float __x) __attribute__ ((__nothrow__ , __leaf__)); extern long int __lrintf (float __x) __attribute__ ((__nothrow__ , __leaf__));
extern long long int llrintf (float __x) __attribute__ ((__nothrow__ , __leaf__)); extern long long int __llrintf (float __x) __attribute__ ((__nothrow__ , __leaf__));


 
extern long int lroundf (float __x) __attribute__ ((__nothrow__ , __leaf__)); extern long int __lroundf (float __x) __attribute__ ((__nothrow__ , __leaf__));
extern long long int llroundf (float __x) __attribute__ ((__nothrow__ , __leaf__)); extern long long int __llroundf (float __x) __attribute__ ((__nothrow__ , __leaf__));


 
extern float fdimf (float __x, float __y) __attribute__ ((__nothrow__ , __leaf__)); extern float __fdimf (float __x, float __y) __attribute__ ((__nothrow__ , __leaf__));

 
extern float fmaxf (float __x, float __y) __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__const__)); extern float __fmaxf (float __x, float __y) __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__const__));

 
extern float fminf (float __x, float __y) __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__const__)); extern float __fminf (float __x, float __y) __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__const__));


 
extern int __fpclassifyf (float __value) __attribute__ ((__nothrow__ , __leaf__))
     __attribute__ ((__const__));

 
extern int __signbitf (float __value) __attribute__ ((__nothrow__ , __leaf__))
     __attribute__ ((__const__));


 
extern float fmaf (float __x, float __y, float __z) __attribute__ ((__nothrow__ , __leaf__)); extern float __fmaf (float __x, float __y, float __z) __attribute__ ((__nothrow__ , __leaf__));



 
extern float scalbf (float __x, float __n) __attribute__ ((__nothrow__ , __leaf__)); extern float __scalbf (float __x, float __n) __attribute__ ((__nothrow__ , __leaf__));



 

















 

























 



 


 
extern long double acosl (long double __x) __attribute__ ((__nothrow__ , __leaf__)); extern long double __acosl (long double __x) __attribute__ ((__nothrow__ , __leaf__));
 
extern long double asinl (long double __x) __attribute__ ((__nothrow__ , __leaf__)); extern long double __asinl (long double __x) __attribute__ ((__nothrow__ , __leaf__));
 
extern long double atanl (long double __x) __attribute__ ((__nothrow__ , __leaf__)); extern long double __atanl (long double __x) __attribute__ ((__nothrow__ , __leaf__));
 
extern long double atan2l (long double __y, long double __x) __attribute__ ((__nothrow__ , __leaf__)); extern long double __atan2l (long double __y, long double __x) __attribute__ ((__nothrow__ , __leaf__));

 
extern long double cosl (long double __x) __attribute__ ((__nothrow__ , __leaf__)); extern long double __cosl (long double __x) __attribute__ ((__nothrow__ , __leaf__));
 
extern long double sinl (long double __x) __attribute__ ((__nothrow__ , __leaf__)); extern long double __sinl (long double __x) __attribute__ ((__nothrow__ , __leaf__));
 
extern long double tanl (long double __x) __attribute__ ((__nothrow__ , __leaf__)); extern long double __tanl (long double __x) __attribute__ ((__nothrow__ , __leaf__));

 

 
extern long double coshl (long double __x) __attribute__ ((__nothrow__ , __leaf__)); extern long double __coshl (long double __x) __attribute__ ((__nothrow__ , __leaf__));
 
extern long double sinhl (long double __x) __attribute__ ((__nothrow__ , __leaf__)); extern long double __sinhl (long double __x) __attribute__ ((__nothrow__ , __leaf__));
 
extern long double tanhl (long double __x) __attribute__ ((__nothrow__ , __leaf__)); extern long double __tanhl (long double __x) __attribute__ ((__nothrow__ , __leaf__));




 
extern long double acoshl (long double __x) __attribute__ ((__nothrow__ , __leaf__)); extern long double __acoshl (long double __x) __attribute__ ((__nothrow__ , __leaf__));
 
extern long double asinhl (long double __x) __attribute__ ((__nothrow__ , __leaf__)); extern long double __asinhl (long double __x) __attribute__ ((__nothrow__ , __leaf__));
 
extern long double atanhl (long double __x) __attribute__ ((__nothrow__ , __leaf__)); extern long double __atanhl (long double __x) __attribute__ ((__nothrow__ , __leaf__));


 


 
extern long double expl (long double __x) __attribute__ ((__nothrow__ , __leaf__)); extern long double __expl (long double __x) __attribute__ ((__nothrow__ , __leaf__));

 
extern long double frexpl (long double __x, int *__exponent) __attribute__ ((__nothrow__ , __leaf__)); extern long double __frexpl (long double __x, int *__exponent) __attribute__ ((__nothrow__ , __leaf__));

 
extern long double ldexpl (long double __x, int __exponent) __attribute__ ((__nothrow__ , __leaf__)); extern long double __ldexpl (long double __x, int __exponent) __attribute__ ((__nothrow__ , __leaf__));

 
extern long double logl (long double __x) __attribute__ ((__nothrow__ , __leaf__)); extern long double __logl (long double __x) __attribute__ ((__nothrow__ , __leaf__));

 
extern long double log10l (long double __x) __attribute__ ((__nothrow__ , __leaf__)); extern long double __log10l (long double __x) __attribute__ ((__nothrow__ , __leaf__));

 
extern long double modfl (long double __x, long double *__iptr) __attribute__ ((__nothrow__ , __leaf__)); extern long double __modfl (long double __x, long double *__iptr) __attribute__ ((__nothrow__ , __leaf__))
     __attribute__ ((__nonnull__ (2)));




 
extern long double expm1l (long double __x) __attribute__ ((__nothrow__ , __leaf__)); extern long double __expm1l (long double __x) __attribute__ ((__nothrow__ , __leaf__));

 
extern long double log1pl (long double __x) __attribute__ ((__nothrow__ , __leaf__)); extern long double __log1pl (long double __x) __attribute__ ((__nothrow__ , __leaf__));

 
extern long double logbl (long double __x) __attribute__ ((__nothrow__ , __leaf__)); extern long double __logbl (long double __x) __attribute__ ((__nothrow__ , __leaf__));



 
extern long double exp2l (long double __x) __attribute__ ((__nothrow__ , __leaf__)); extern long double __exp2l (long double __x) __attribute__ ((__nothrow__ , __leaf__));

 
extern long double log2l (long double __x) __attribute__ ((__nothrow__ , __leaf__)); extern long double __log2l (long double __x) __attribute__ ((__nothrow__ , __leaf__));



 


 
extern long double powl (long double __x, long double __y) __attribute__ ((__nothrow__ , __leaf__)); extern long double __powl (long double __x, long double __y) __attribute__ ((__nothrow__ , __leaf__));

 
extern long double sqrtl (long double __x) __attribute__ ((__nothrow__ , __leaf__)); extern long double __sqrtl (long double __x) __attribute__ ((__nothrow__ , __leaf__));



 
extern long double hypotl (long double __x, long double __y) __attribute__ ((__nothrow__ , __leaf__)); extern long double __hypotl (long double __x, long double __y) __attribute__ ((__nothrow__ , __leaf__));



 
extern long double cbrtl (long double __x) __attribute__ ((__nothrow__ , __leaf__)); extern long double __cbrtl (long double __x) __attribute__ ((__nothrow__ , __leaf__));



 


 
extern long double ceill (long double __x) __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__const__)); extern long double __ceill (long double __x) __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__const__));

 
extern long double fabsl (long double __x) __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__const__)); extern long double __fabsl (long double __x) __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__const__));

 
extern long double floorl (long double __x) __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__const__)); extern long double __floorl (long double __x) __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__const__));

 
extern long double fmodl (long double __x, long double __y) __attribute__ ((__nothrow__ , __leaf__)); extern long double __fmodl (long double __x, long double __y) __attribute__ ((__nothrow__ , __leaf__));



 
extern int __isinfl (long double __value) __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__const__));

 
extern int __finitel (long double __value) __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__const__));



 
extern int isinfl (long double __value) __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__const__));

 
extern int finitel (long double __value) __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__const__));

 
extern long double dreml (long double __x, long double __y) __attribute__ ((__nothrow__ , __leaf__)); extern long double __dreml (long double __x, long double __y) __attribute__ ((__nothrow__ , __leaf__));


 
extern long double significandl (long double __x) __attribute__ ((__nothrow__ , __leaf__)); extern long double __significandl (long double __x) __attribute__ ((__nothrow__ , __leaf__));


 
extern long double copysignl (long double __x, long double __y) __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__const__)); extern long double __copysignl (long double __x, long double __y) __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__const__));



 
extern long double nanl (const char *__tagb) __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__const__)); extern long double __nanl (const char *__tagb) __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__const__));



 
extern int __isnanl (long double __value) __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__const__));

 
extern int isnanl (long double __value) __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__const__));

 
extern long double j0l (long double) __attribute__ ((__nothrow__ , __leaf__)); extern long double __j0l (long double) __attribute__ ((__nothrow__ , __leaf__));
extern long double j1l (long double) __attribute__ ((__nothrow__ , __leaf__)); extern long double __j1l (long double) __attribute__ ((__nothrow__ , __leaf__));
extern long double jnl (int, long double) __attribute__ ((__nothrow__ , __leaf__)); extern long double __jnl (int, long double) __attribute__ ((__nothrow__ , __leaf__));
extern long double y0l (long double) __attribute__ ((__nothrow__ , __leaf__)); extern long double __y0l (long double) __attribute__ ((__nothrow__ , __leaf__));
extern long double y1l (long double) __attribute__ ((__nothrow__ , __leaf__)); extern long double __y1l (long double) __attribute__ ((__nothrow__ , __leaf__));
extern long double ynl (int, long double) __attribute__ ((__nothrow__ , __leaf__)); extern long double __ynl (int, long double) __attribute__ ((__nothrow__ , __leaf__));



 
extern long double erfl (long double) __attribute__ ((__nothrow__ , __leaf__)); extern long double __erfl (long double) __attribute__ ((__nothrow__ , __leaf__));
extern long double erfcl (long double) __attribute__ ((__nothrow__ , __leaf__)); extern long double __erfcl (long double) __attribute__ ((__nothrow__ , __leaf__));
extern long double lgammal (long double) __attribute__ ((__nothrow__ , __leaf__)); extern long double __lgammal (long double) __attribute__ ((__nothrow__ , __leaf__));



 
extern long double tgammal (long double) __attribute__ ((__nothrow__ , __leaf__)); extern long double __tgammal (long double) __attribute__ ((__nothrow__ , __leaf__));


 
extern long double gammal (long double) __attribute__ ((__nothrow__ , __leaf__)); extern long double __gammal (long double) __attribute__ ((__nothrow__ , __leaf__));



 
extern long double lgammal_r (long double, int *__signgamp) __attribute__ ((__nothrow__ , __leaf__)); extern long double __lgammal_r (long double, int *__signgamp) __attribute__ ((__nothrow__ , __leaf__));




 
extern long double rintl (long double __x) __attribute__ ((__nothrow__ , __leaf__)); extern long double __rintl (long double __x) __attribute__ ((__nothrow__ , __leaf__));

 
extern long double nextafterl (long double __x, long double __y) __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__const__)); extern long double __nextafterl (long double __x, long double __y) __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__const__));
extern long double nexttowardl (long double __x, long double __y) __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__const__)); extern long double __nexttowardl (long double __x, long double __y) __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__const__));

 
extern long double remainderl (long double __x, long double __y) __attribute__ ((__nothrow__ , __leaf__)); extern long double __remainderl (long double __x, long double __y) __attribute__ ((__nothrow__ , __leaf__));

 
extern long double scalbnl (long double __x, int __n) __attribute__ ((__nothrow__ , __leaf__)); extern long double __scalbnl (long double __x, int __n) __attribute__ ((__nothrow__ , __leaf__));

 
extern int ilogbl (long double __x) __attribute__ ((__nothrow__ , __leaf__)); extern int __ilogbl (long double __x) __attribute__ ((__nothrow__ , __leaf__));

 
extern long double scalblnl (long double __x, long int __n) __attribute__ ((__nothrow__ , __leaf__)); extern long double __scalblnl (long double __x, long int __n) __attribute__ ((__nothrow__ , __leaf__));


 
extern long double nearbyintl (long double __x) __attribute__ ((__nothrow__ , __leaf__)); extern long double __nearbyintl (long double __x) __attribute__ ((__nothrow__ , __leaf__));


 
extern long double roundl (long double __x) __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__const__)); extern long double __roundl (long double __x) __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__const__));


 
extern long double truncl (long double __x) __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__const__)); extern long double __truncl (long double __x) __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__const__));



 
extern long double remquol (long double __x, long double __y, int *__quo) __attribute__ ((__nothrow__ , __leaf__)); extern long double __remquol (long double __x, long double __y, int *__quo) __attribute__ ((__nothrow__ , __leaf__));


 


 
extern long int lrintl (long double __x) __attribute__ ((__nothrow__ , __leaf__)); extern long int __lrintl (long double __x) __attribute__ ((__nothrow__ , __leaf__));
extern long long int llrintl (long double __x) __attribute__ ((__nothrow__ , __leaf__)); extern long long int __llrintl (long double __x) __attribute__ ((__nothrow__ , __leaf__));


 
extern long int lroundl (long double __x) __attribute__ ((__nothrow__ , __leaf__)); extern long int __lroundl (long double __x) __attribute__ ((__nothrow__ , __leaf__));
extern long long int llroundl (long double __x) __attribute__ ((__nothrow__ , __leaf__)); extern long long int __llroundl (long double __x) __attribute__ ((__nothrow__ , __leaf__));


 
extern long double fdiml (long double __x, long double __y) __attribute__ ((__nothrow__ , __leaf__)); extern long double __fdiml (long double __x, long double __y) __attribute__ ((__nothrow__ , __leaf__));

 
extern long double fmaxl (long double __x, long double __y) __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__const__)); extern long double __fmaxl (long double __x, long double __y) __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__const__));

 
extern long double fminl (long double __x, long double __y) __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__const__)); extern long double __fminl (long double __x, long double __y) __attribute__ ((__nothrow__ , __leaf__)) __attribute__ ((__const__));


 
extern int __fpclassifyl (long double __value) __attribute__ ((__nothrow__ , __leaf__))
     __attribute__ ((__const__));

 
extern int __signbitl (long double __value) __attribute__ ((__nothrow__ , __leaf__))
     __attribute__ ((__const__));


 
extern long double fmal (long double __x, long double __y, long double __z) __attribute__ ((__nothrow__ , __leaf__)); extern long double __fmal (long double __x, long double __y, long double __z) __attribute__ ((__nothrow__ , __leaf__));



 
extern long double scalbl (long double __x, long double __n) __attribute__ ((__nothrow__ , __leaf__)); extern long double __scalbl (long double __x, long double __n) __attribute__ ((__nothrow__ , __leaf__));




 
extern int signgam;


 
































 

 
enum
  {
    FP_NAN =
      0,
    FP_INFINITE =
      1,
    FP_ZERO =
      2,
    FP_SUBNORMAL =
      3,
    FP_NORMAL =
      4
  };

 

 

 

 


 

 

 



 


 
typedef enum
{
  _IEEE_ = -1,	 
  _SVID_,	 
  _XOPEN_,	 
  _POSIX_,
  _ISOC_	 
} _LIB_VERSION_TYPE;



 
extern _LIB_VERSION_TYPE _LIB_VERSION;






 
struct exception
  {
    int type;
    char *name;
    double arg1;
    double arg2;
    double retval;
  };

extern int matherr (struct exception *__exc);


 

 



 



 




 






 

 
















 




 


 


 
extern __inline __attribute__ ((__always_inline__)) __attribute__ ((__gnu_inline__)) int
__attribute__ ((__nothrow__ , __leaf__)) __signbitf (float __x)
{
  int __m;
  __asm ("pmovmskb %1, %0" : "=r" (__m) : "x" (__x));
  return (__m & 0x8) != 0;
}
extern __inline __attribute__ ((__always_inline__)) __attribute__ ((__gnu_inline__)) int
__attribute__ ((__nothrow__ , __leaf__)) __signbit (double __x)
{
  int __m;
  __asm ("pmovmskb %1, %0" : "=r" (__m) : "x" (__x));
  return (__m & 0x80) != 0;
}
extern __inline __attribute__ ((__always_inline__)) __attribute__ ((__gnu_inline__)) int
__attribute__ ((__nothrow__ , __leaf__)) __signbitl (long double __x)
{
  __extension__ union { long double __l; int __i[3]; } __u = { __l: __x };
  return (__u.__i[2] & 0x8000) != 0;
}





 



 

 

 

 

 

 

 

 






 


 








 











 




 



 extern int fpclassifyf    ( float              __x ) ;
 extern int fpclassify     ( double           __x ) ;
 extern int fpclassifyd    ( double             __x ) ;
 extern int fpclassifyl    ( long double            __x ) ;

 extern int isinff         ( float              __x ) ;
 extern int isinf          ( double           __x ) ;
 extern int isinfd         ( double             __x ) ;
 extern int isinfl         ( long double            __x ) ;

 extern int isnanf         ( float              __x ) ;
 extern int isnan          ( double           __x ) ;
 extern int isnand         ( double             __x ) ;
 extern int isnanl         ( long double            __x ) ;

 extern int isnormalf      ( float              __x ) ;
 extern int isnormal       ( double           __x ) ;
 extern int isnormald      ( double             __x ) ;
 extern int isnormall      ( long double            __x ) ;

 extern int isfinitef      ( float              __x ) ;
 extern int isfinite       ( double           __x ) ;
 extern int isfinited      ( double             __x ) ;
 extern int isfinitel      ( long double            __x ) ;

 extern int finitef        ( float              __x ) ;
 extern int finite         ( double             __x ) ;
 extern int finited        ( double             __x ) ;
 extern int finitel        ( long double            __x ) ;

 extern int signbitf       ( float              __x ) ;
 extern int signbit        ( double             __x ) ;
 extern int signbitd       ( double             __x ) ;
 extern int signbitl       ( long double            __x ) ;


 extern int __fpclassifyf  ( float              __x ) ;
 extern int __fpclassify   ( double           __x ) ;
 extern int __fpclassifyd  ( double             __x ) ;
 extern int __fpclassifyl  ( long double            __x ) ;

 extern int __isinff       ( float              __x ) ;
 extern int __isinf        ( double           __x ) ;
 extern int __isinfd       ( double             __x ) ;
 extern int __isinfl       ( long double            __x ) ;

 extern int __isnanf       ( float              __x ) ;
 extern int __isnan        ( double           __x ) ;
 extern int __isnand       ( double             __x ) ;
 extern int __isnanl       ( long double            __x ) ;

 extern int __isnormalf    ( float              __x ) ;
 extern int __isnormal     ( double           __x ) ;
 extern int __isnormald    ( double             __x ) ;
 extern int __isnormall    ( long double            __x ) ;

 extern int __isfinitef    ( float              __x ) ;
 extern int __isfinite     ( double           __x ) ;
 extern int __isfinited    ( double             __x ) ;
 extern int __isfinitel    ( long double            __x ) ;

 extern int __finitef      ( float              __x ) ;
 extern int __finite       ( double             __x ) ;
 extern int __finited      ( double             __x ) ;
 extern int __finitel      ( long double            __x ) ;

 extern int __signbitf     ( float              __x ) ;
 extern int __signbit      ( double             __x ) ;
 extern int __signbitd     ( double             __x ) ;
 extern int __signbitl     ( long double            __x ) ;









 

 extern int isgreaterf( float __xf, float __yf );
 extern int isgreater( double __xd, double __yd );
 extern int isgreaterl( long double __xl, long double __yl );
 extern int __isgreaterf( float __xf, float __yf );
 extern int __isgreater( double __xd, double __yd );
 extern int __isgreaterl( long double __xl, long double __yl );

 extern int isgreaterequalf( float __xf, float __yf );
 extern int isgreaterequal( double __xd, double __yd );
 extern int isgreaterequall( long double __xl, long double __yl );
 extern int __isgreaterequalf( float __xf, float __yf );
 extern int __isgreaterequal( double __xd, double __yd );
 extern int __isgreaterequall( long double __xl, long double __yl );

 extern int islessf( float __xf, float __yf );
 extern int isless( double __xd, double __yd );
 extern int islessl( long double __xl, long double __yl );
 extern int __islessf( float __xf, float __yf );
 extern int __isless( double __xd, double __yd );
 extern int __islessl( long double __xl, long double __yl );

 int islessequalf( float __xf, float __yf );
 extern int islessequal( double __xd, double __yd );
 extern int islessequall( long double __xl, long double __yl );
 extern int __islessequalf( float __xf, float __yf );
 extern int __islessequal( double __xd, double __yd );
 extern int __islessequall( long double __xl, long double __yl );

 extern int islessgreaterf( float __xf, float __yf );
 extern int islessgreater( double __xd, double __yd );
 extern int islessgreaterl( long double __xl, long double __yl );
 extern int __islessgreaterf( float __xf, float __yf );
 extern int __islessgreater( double __xd, double __yd );
 extern int __islessgreaterl( long double __xl, long double __yl );

 extern int isunorderedf( float __xf, float __yf );
 extern int isunordered( double __xd, double __yd );
 extern int isunorderedl( long double __xl, long double __yl );
 extern int __isunorderedf( float __xf, float __yf );
 extern int __isunordered( double __xd, double __yd );
 extern int __isunorderedl( long double __xl, long double __yl );





 

 




 




 




 
 
 
 
 








 







 





 


  
extern double    gamma( double __x );
extern float     gammaf( float __x );

 
extern double    lgamma_r(double __x, int *__signgam);
extern float     lgammaf_r( float __x, int *__signgam );

extern double    gamma_r( double __x, int *__signgam );
extern float     gammaf_r( float __x, int *__signgam );


 







 





 





 








 

 


typedef struct ____exception {
    int     type;
    const char  *name;
    double  arg1;
    double  arg2;
    double  retval;
} ___exception;


typedef struct ____exceptionf {
    int    type;
    const char *name;
    float  arg1;
    float  arg2;
    float  retval;
} ___exceptionf;

typedef struct ____exceptionl {
    int      type;
    const char   *name;
    long double  arg1;
    long double  arg2;
    long double  retval;
} ___exceptionl;

 extern int  matherrf( struct ____exceptionf *__e );
 extern int  matherrl( struct ____exceptionl *__e );
















 

typedef int (  *___pmatherr )( struct ____exception  *__e );
typedef int (  *___pmatherrf )( struct ____exceptionf *__e );
typedef int (  *___pmatherrl )( struct ____exceptionl *__e );

 extern ___pmatherr   __libm_setusermatherr( ___pmatherr  __user_matherr );
 extern ___pmatherrf  __libm_setusermatherrf( ___pmatherrf __user_matherrf );
 extern ___pmatherrl  __libm_setusermatherrl( ___pmatherrl __user_matherrl );

 


extern _LIB_VERSION_TYPE  _LIB_VERSIONIMF;


 








 




 








 


 
















 



 



 















 







 




 


 

 




 



















 






















 




 

 
extern double _Complex cacos (double _Complex __z) __attribute__ ((__nothrow__ , __leaf__)); extern double _Complex __cacos (double _Complex __z) __attribute__ ((__nothrow__ , __leaf__));
 
extern double _Complex casin (double _Complex __z) __attribute__ ((__nothrow__ , __leaf__)); extern double _Complex __casin (double _Complex __z) __attribute__ ((__nothrow__ , __leaf__));
 
extern double _Complex catan (double _Complex __z) __attribute__ ((__nothrow__ , __leaf__)); extern double _Complex __catan (double _Complex __z) __attribute__ ((__nothrow__ , __leaf__));

 
extern double _Complex ccos (double _Complex __z) __attribute__ ((__nothrow__ , __leaf__)); extern double _Complex __ccos (double _Complex __z) __attribute__ ((__nothrow__ , __leaf__));
 
extern double _Complex csin (double _Complex __z) __attribute__ ((__nothrow__ , __leaf__)); extern double _Complex __csin (double _Complex __z) __attribute__ ((__nothrow__ , __leaf__));
 
extern double _Complex ctan (double _Complex __z) __attribute__ ((__nothrow__ , __leaf__)); extern double _Complex __ctan (double _Complex __z) __attribute__ ((__nothrow__ , __leaf__));


 

 
extern double _Complex cacosh (double _Complex __z) __attribute__ ((__nothrow__ , __leaf__)); extern double _Complex __cacosh (double _Complex __z) __attribute__ ((__nothrow__ , __leaf__));
 
extern double _Complex casinh (double _Complex __z) __attribute__ ((__nothrow__ , __leaf__)); extern double _Complex __casinh (double _Complex __z) __attribute__ ((__nothrow__ , __leaf__));
 
extern double _Complex catanh (double _Complex __z) __attribute__ ((__nothrow__ , __leaf__)); extern double _Complex __catanh (double _Complex __z) __attribute__ ((__nothrow__ , __leaf__));

 
extern double _Complex ccosh (double _Complex __z) __attribute__ ((__nothrow__ , __leaf__)); extern double _Complex __ccosh (double _Complex __z) __attribute__ ((__nothrow__ , __leaf__));
 
extern double _Complex csinh (double _Complex __z) __attribute__ ((__nothrow__ , __leaf__)); extern double _Complex __csinh (double _Complex __z) __attribute__ ((__nothrow__ , __leaf__));
 
extern double _Complex ctanh (double _Complex __z) __attribute__ ((__nothrow__ , __leaf__)); extern double _Complex __ctanh (double _Complex __z) __attribute__ ((__nothrow__ , __leaf__));


 

 
extern double _Complex cexp (double _Complex __z) __attribute__ ((__nothrow__ , __leaf__)); extern double _Complex __cexp (double _Complex __z) __attribute__ ((__nothrow__ , __leaf__));

 
extern double _Complex clog (double _Complex __z) __attribute__ ((__nothrow__ , __leaf__)); extern double _Complex __clog (double _Complex __z) __attribute__ ((__nothrow__ , __leaf__));


 

 
extern double _Complex cpow (double _Complex __x, double _Complex __y) __attribute__ ((__nothrow__ , __leaf__)); extern double _Complex __cpow (double _Complex __x, double _Complex __y) __attribute__ ((__nothrow__ , __leaf__));

 
extern double _Complex csqrt (double _Complex __z) __attribute__ ((__nothrow__ , __leaf__)); extern double _Complex __csqrt (double _Complex __z) __attribute__ ((__nothrow__ , __leaf__));


 

 
extern double cabs (double _Complex __z) __attribute__ ((__nothrow__ , __leaf__)); extern double __cabs (double _Complex __z) __attribute__ ((__nothrow__ , __leaf__));

 
extern double carg (double _Complex __z) __attribute__ ((__nothrow__ , __leaf__)); extern double __carg (double _Complex __z) __attribute__ ((__nothrow__ , __leaf__));

 
extern double _Complex conj (double _Complex __z) __attribute__ ((__nothrow__ , __leaf__)); extern double _Complex __conj (double _Complex __z) __attribute__ ((__nothrow__ , __leaf__));

 
extern double _Complex cproj (double _Complex __z) __attribute__ ((__nothrow__ , __leaf__)); extern double _Complex __cproj (double _Complex __z) __attribute__ ((__nothrow__ , __leaf__));


 

 
extern double cimag (double _Complex __z) __attribute__ ((__nothrow__ , __leaf__)); extern double __cimag (double _Complex __z) __attribute__ ((__nothrow__ , __leaf__));

 
extern double creal (double _Complex __z) __attribute__ ((__nothrow__ , __leaf__)); extern double __creal (double _Complex __z) __attribute__ ((__nothrow__ , __leaf__));




 

 

















 






















 




 

 
extern float _Complex cacosf (float _Complex __z) __attribute__ ((__nothrow__ , __leaf__)); extern float _Complex __cacosf (float _Complex __z) __attribute__ ((__nothrow__ , __leaf__));
 
extern float _Complex casinf (float _Complex __z) __attribute__ ((__nothrow__ , __leaf__)); extern float _Complex __casinf (float _Complex __z) __attribute__ ((__nothrow__ , __leaf__));
 
extern float _Complex catanf (float _Complex __z) __attribute__ ((__nothrow__ , __leaf__)); extern float _Complex __catanf (float _Complex __z) __attribute__ ((__nothrow__ , __leaf__));

 
extern float _Complex ccosf (float _Complex __z) __attribute__ ((__nothrow__ , __leaf__)); extern float _Complex __ccosf (float _Complex __z) __attribute__ ((__nothrow__ , __leaf__));
 
extern float _Complex csinf (float _Complex __z) __attribute__ ((__nothrow__ , __leaf__)); extern float _Complex __csinf (float _Complex __z) __attribute__ ((__nothrow__ , __leaf__));
 
extern float _Complex ctanf (float _Complex __z) __attribute__ ((__nothrow__ , __leaf__)); extern float _Complex __ctanf (float _Complex __z) __attribute__ ((__nothrow__ , __leaf__));


 

 
extern float _Complex cacoshf (float _Complex __z) __attribute__ ((__nothrow__ , __leaf__)); extern float _Complex __cacoshf (float _Complex __z) __attribute__ ((__nothrow__ , __leaf__));
 
extern float _Complex casinhf (float _Complex __z) __attribute__ ((__nothrow__ , __leaf__)); extern float _Complex __casinhf (float _Complex __z) __attribute__ ((__nothrow__ , __leaf__));
 
extern float _Complex catanhf (float _Complex __z) __attribute__ ((__nothrow__ , __leaf__)); extern float _Complex __catanhf (float _Complex __z) __attribute__ ((__nothrow__ , __leaf__));

 
extern float _Complex ccoshf (float _Complex __z) __attribute__ ((__nothrow__ , __leaf__)); extern float _Complex __ccoshf (float _Complex __z) __attribute__ ((__nothrow__ , __leaf__));
 
extern float _Complex csinhf (float _Complex __z) __attribute__ ((__nothrow__ , __leaf__)); extern float _Complex __csinhf (float _Complex __z) __attribute__ ((__nothrow__ , __leaf__));
 
extern float _Complex ctanhf (float _Complex __z) __attribute__ ((__nothrow__ , __leaf__)); extern float _Complex __ctanhf (float _Complex __z) __attribute__ ((__nothrow__ , __leaf__));


 

 
extern float _Complex cexpf (float _Complex __z) __attribute__ ((__nothrow__ , __leaf__)); extern float _Complex __cexpf (float _Complex __z) __attribute__ ((__nothrow__ , __leaf__));

 
extern float _Complex clogf (float _Complex __z) __attribute__ ((__nothrow__ , __leaf__)); extern float _Complex __clogf (float _Complex __z) __attribute__ ((__nothrow__ , __leaf__));


 

 
extern float _Complex cpowf (float _Complex __x, float _Complex __y) __attribute__ ((__nothrow__ , __leaf__)); extern float _Complex __cpowf (float _Complex __x, float _Complex __y) __attribute__ ((__nothrow__ , __leaf__));

 
extern float _Complex csqrtf (float _Complex __z) __attribute__ ((__nothrow__ , __leaf__)); extern float _Complex __csqrtf (float _Complex __z) __attribute__ ((__nothrow__ , __leaf__));


 

 
extern float cabsf (float _Complex __z) __attribute__ ((__nothrow__ , __leaf__)); extern float __cabsf (float _Complex __z) __attribute__ ((__nothrow__ , __leaf__));

 
extern float cargf (float _Complex __z) __attribute__ ((__nothrow__ , __leaf__)); extern float __cargf (float _Complex __z) __attribute__ ((__nothrow__ , __leaf__));

 
extern float _Complex conjf (float _Complex __z) __attribute__ ((__nothrow__ , __leaf__)); extern float _Complex __conjf (float _Complex __z) __attribute__ ((__nothrow__ , __leaf__));

 
extern float _Complex cprojf (float _Complex __z) __attribute__ ((__nothrow__ , __leaf__)); extern float _Complex __cprojf (float _Complex __z) __attribute__ ((__nothrow__ , __leaf__));


 

 
extern float cimagf (float _Complex __z) __attribute__ ((__nothrow__ , __leaf__)); extern float __cimagf (float _Complex __z) __attribute__ ((__nothrow__ , __leaf__));

 
extern float crealf (float _Complex __z) __attribute__ ((__nothrow__ , __leaf__)); extern float __crealf (float _Complex __z) __attribute__ ((__nothrow__ , __leaf__));




 


 


















 






















 




 

 
extern long double _Complex cacosl (long double _Complex __z) __attribute__ ((__nothrow__ , __leaf__)); extern long double _Complex __cacosl (long double _Complex __z) __attribute__ ((__nothrow__ , __leaf__));
 
extern long double _Complex casinl (long double _Complex __z) __attribute__ ((__nothrow__ , __leaf__)); extern long double _Complex __casinl (long double _Complex __z) __attribute__ ((__nothrow__ , __leaf__));
 
extern long double _Complex catanl (long double _Complex __z) __attribute__ ((__nothrow__ , __leaf__)); extern long double _Complex __catanl (long double _Complex __z) __attribute__ ((__nothrow__ , __leaf__));

 
extern long double _Complex ccosl (long double _Complex __z) __attribute__ ((__nothrow__ , __leaf__)); extern long double _Complex __ccosl (long double _Complex __z) __attribute__ ((__nothrow__ , __leaf__));
 
extern long double _Complex csinl (long double _Complex __z) __attribute__ ((__nothrow__ , __leaf__)); extern long double _Complex __csinl (long double _Complex __z) __attribute__ ((__nothrow__ , __leaf__));
 
extern long double _Complex ctanl (long double _Complex __z) __attribute__ ((__nothrow__ , __leaf__)); extern long double _Complex __ctanl (long double _Complex __z) __attribute__ ((__nothrow__ , __leaf__));


 

 
extern long double _Complex cacoshl (long double _Complex __z) __attribute__ ((__nothrow__ , __leaf__)); extern long double _Complex __cacoshl (long double _Complex __z) __attribute__ ((__nothrow__ , __leaf__));
 
extern long double _Complex casinhl (long double _Complex __z) __attribute__ ((__nothrow__ , __leaf__)); extern long double _Complex __casinhl (long double _Complex __z) __attribute__ ((__nothrow__ , __leaf__));
 
extern long double _Complex catanhl (long double _Complex __z) __attribute__ ((__nothrow__ , __leaf__)); extern long double _Complex __catanhl (long double _Complex __z) __attribute__ ((__nothrow__ , __leaf__));

 
extern long double _Complex ccoshl (long double _Complex __z) __attribute__ ((__nothrow__ , __leaf__)); extern long double _Complex __ccoshl (long double _Complex __z) __attribute__ ((__nothrow__ , __leaf__));
 
extern long double _Complex csinhl (long double _Complex __z) __attribute__ ((__nothrow__ , __leaf__)); extern long double _Complex __csinhl (long double _Complex __z) __attribute__ ((__nothrow__ , __leaf__));
 
extern long double _Complex ctanhl (long double _Complex __z) __attribute__ ((__nothrow__ , __leaf__)); extern long double _Complex __ctanhl (long double _Complex __z) __attribute__ ((__nothrow__ , __leaf__));


 

 
extern long double _Complex cexpl (long double _Complex __z) __attribute__ ((__nothrow__ , __leaf__)); extern long double _Complex __cexpl (long double _Complex __z) __attribute__ ((__nothrow__ , __leaf__));

 
extern long double _Complex clogl (long double _Complex __z) __attribute__ ((__nothrow__ , __leaf__)); extern long double _Complex __clogl (long double _Complex __z) __attribute__ ((__nothrow__ , __leaf__));


 

 
extern long double _Complex cpowl (long double _Complex __x, long double _Complex __y) __attribute__ ((__nothrow__ , __leaf__)); extern long double _Complex __cpowl (long double _Complex __x, long double _Complex __y) __attribute__ ((__nothrow__ , __leaf__));

 
extern long double _Complex csqrtl (long double _Complex __z) __attribute__ ((__nothrow__ , __leaf__)); extern long double _Complex __csqrtl (long double _Complex __z) __attribute__ ((__nothrow__ , __leaf__));


 

 
extern long double cabsl (long double _Complex __z) __attribute__ ((__nothrow__ , __leaf__)); extern long double __cabsl (long double _Complex __z) __attribute__ ((__nothrow__ , __leaf__));

 
extern long double cargl (long double _Complex __z) __attribute__ ((__nothrow__ , __leaf__)); extern long double __cargl (long double _Complex __z) __attribute__ ((__nothrow__ , __leaf__));

 
extern long double _Complex conjl (long double _Complex __z) __attribute__ ((__nothrow__ , __leaf__)); extern long double _Complex __conjl (long double _Complex __z) __attribute__ ((__nothrow__ , __leaf__));

 
extern long double _Complex cprojl (long double _Complex __z) __attribute__ ((__nothrow__ , __leaf__)); extern long double _Complex __cprojl (long double _Complex __z) __attribute__ ((__nothrow__ , __leaf__));


 

 
extern long double cimagl (long double _Complex __z) __attribute__ ((__nothrow__ , __leaf__)); extern long double __cimagl (long double _Complex __z) __attribute__ ((__nothrow__ , __leaf__));

 
extern long double creall (long double _Complex __z) __attribute__ ((__nothrow__ , __leaf__)); extern long double __creall (long double _Complex __z) __attribute__ ((__nothrow__ , __leaf__));




 





 


 








 







 

 



 


extern double _Complex  cis( double __x );
extern float _Complex  cisf( float __x );
extern long double _Complex  cisl( long double __x );

extern double _Complex  cisd( double __x );
extern float _Complex  cisdf( float __x );
extern long double _Complex  cisdl( long double __x );

 


extern double _Complex  cexp2( double _Complex __z );
extern float _Complex  cexp2f( float _Complex __z );
extern long double _Complex  cexp2l( long double _Complex __z );

extern double _Complex  cexp10( double _Complex __z );
extern float _Complex  cexp10f( float _Complex __z );
extern long double _Complex  cexp10l( long double _Complex __z );


 


extern double _Complex  clog2( double _Complex __z );
extern float _Complex  clog2f( float _Complex __z );
extern long double _Complex  clog2l( long double _Complex __z );

extern double _Complex  clog10( double _Complex __z );
extern float _Complex  clog10f( float _Complex __z );
extern long double _Complex  clog10l( long double _Complex __z );

 


 






 








 








int numroc_(int*, int*, int*, int*, int*);


void pdlacpy_(char*, int*, int*, double*, int*, int*, int*, double*, int*, int*, int*);
void pdtran_(int*, int*, double*, double*, int*, int*, int*, double*, double*, int*, int*, int*);


void pslacpy_(char*, int*, int*, float*, int*, int*, int*, float*, int*, int*, int*);
void pstran_(int*, int*, float*, float*, int*, int*, int*, float*, float*, int*, int*, int*);



void pzlacpy_(char*, int*, int*, double _Complex*, int*, int*, int*, double _Complex*, int*, int*, int*);
void pztranc_(int*, int*, double _Complex*, double _Complex*, int*, int*, int*, double _Complex*, double _Complex*, int*, int*, int*);


void pclacpy_(char*, int*, int*, float _Complex*, int*, int*, int*, float _Complex*, int*, int*, int*);
void pctranc_(int*, int*, float _Complex*, float _Complex*, int*, int*, int*, float _Complex*, float _Complex*, int*, int*, int*);



















 
 










































 
 

 


 



 

 





























 









 




































 



 


















 


 


 


 


 

 

 

 

 
typedef unsigned char		uint8_t;
typedef unsigned short int	uint16_t;
typedef unsigned int		uint32_t;
typedef unsigned long int	uint64_t;


 

 
typedef signed char		int_least8_t;
typedef short int		int_least16_t;
typedef int			int_least32_t;
typedef long int		int_least64_t;

 
typedef unsigned char		uint_least8_t;
typedef unsigned short int	uint_least16_t;
typedef unsigned int		uint_least32_t;
typedef unsigned long int	uint_least64_t;


 

 
typedef signed char		int_fast8_t;
typedef long int		int_fast16_t;
typedef long int		int_fast32_t;
typedef long int		int_fast64_t;

 
typedef unsigned char		uint_fast8_t;
typedef unsigned long int	uint_fast16_t;
typedef unsigned long int	uint_fast32_t;
typedef unsigned long int	uint_fast64_t;


 
typedef long int		intptr_t;
typedef unsigned long int	uintptr_t;


 
typedef long int		intmax_t;
typedef unsigned long int	uintmax_t;



 

 
 

 


 
 

 


 
 

 


 


 
 

 


 

 

 

 

 
 

 


 

 

 




 




 

 

typedef int MPI_Datatype;










 


 






 

 

 


 
 


 

 

 

 
typedef int MPI_Comm;

 
typedef int MPI_Group;

 
typedef int MPI_Win;

 
 

 
typedef struct ADIOI_FileD *MPI_File;

 
typedef int MPI_Op;


 











 


 


 

 

 
typedef enum MPIR_Win_flavor {
    MPI_WIN_FLAVOR_CREATE      = 1,
    MPI_WIN_FLAVOR_ALLOCATE    = 2,
    MPI_WIN_FLAVOR_DYNAMIC     = 3,
    MPI_WIN_FLAVOR_SHARED      = 4
} MPIR_Win_flavor_t;

 
typedef enum MPIR_Win_model {
    MPI_WIN_SEPARATE   = 1,
    MPI_WIN_UNIFIED    = 2
} MPIR_Win_model_t;

 

 
typedef enum MPIR_Topo_type { MPI_GRAPH=1, MPI_CART=2, MPI_DIST_GRAPH=3 } MPIR_Topo_type;

extern  int * const MPI_UNWEIGHTED ;
extern  int * const MPI_WEIGHTS_EMPTY ;



 
typedef void (MPI_Handler_function) ( MPI_Comm *, int *, ... );
typedef int (MPI_Comm_copy_attr_function)(MPI_Comm, int, void *, void *, 
					  void *, int *);
typedef int (MPI_Comm_delete_attr_function)(MPI_Comm, int, void *, void *);
typedef int (MPI_Type_copy_attr_function)(MPI_Datatype, int, void *, void *, 
					  void *, int *);
typedef int (MPI_Type_delete_attr_function)(MPI_Datatype, int, void *, void *);
typedef int (MPI_Win_copy_attr_function)(MPI_Win, int, void *, void *, void *,
					 int *);
typedef int (MPI_Win_delete_attr_function)(MPI_Win, int, void *, void *);
 
typedef void (MPI_Comm_errhandler_function)(MPI_Comm *, int *, ...);
typedef void (MPI_File_errhandler_function)(MPI_File *, int *, ...);
typedef void (MPI_Win_errhandler_function)(MPI_Win *, int *, ...);
 
typedef MPI_Comm_errhandler_function MPI_Comm_errhandler_fn;
typedef MPI_File_errhandler_function MPI_File_errhandler_fn;
typedef MPI_Win_errhandler_function MPI_Win_errhandler_fn;


 



 
typedef int MPI_Errhandler;



 
 
 

 
typedef int MPI_Request;

 
typedef int MPI_Message;

 
typedef void (MPI_User_function) ( void *, void *, int *, MPI_Datatype * ); 

 
typedef int (MPI_Copy_function) ( MPI_Comm, int, void *, void *, void *, int * );
typedef int (MPI_Delete_function) ( MPI_Comm, int, void *, void * );

























 


























 

 
enum MPIR_Combiner_enum {
    MPI_COMBINER_NAMED            = 1,
    MPI_COMBINER_DUP              = 2,
    MPI_COMBINER_CONTIGUOUS       = 3, 
    MPI_COMBINER_VECTOR           = 4,
    MPI_COMBINER_HVECTOR_INTEGER  = 5,
    MPI_COMBINER_HVECTOR          = 6,
    MPI_COMBINER_INDEXED          = 7,
    MPI_COMBINER_HINDEXED_INTEGER = 8, 
    MPI_COMBINER_HINDEXED         = 9, 
    MPI_COMBINER_INDEXED_BLOCK    = 10, 
    MPI_COMBINER_STRUCT_INTEGER   = 11,
    MPI_COMBINER_STRUCT           = 12,
    MPI_COMBINER_SUBARRAY         = 13,
    MPI_COMBINER_DARRAY           = 14,
    MPI_COMBINER_F90_REAL         = 15,
    MPI_COMBINER_F90_COMPLEX      = 16,
    MPI_COMBINER_F90_INTEGER      = 17,
    MPI_COMBINER_RESIZED          = 18,
    MPI_COMBINER_HINDEXED_BLOCK   = 19
};

 
typedef int MPI_Info;

 


 

 

 

 
typedef long MPI_Aint;
typedef int MPI_Fint;
typedef long long MPI_Count;




 

 


 
typedef long long MPI_Offset;



 
typedef struct MPI_Status {
    int count_lo;
    int count_hi_and_cancelled;
    int MPI_SOURCE;
    int MPI_TAG;
    int MPI_ERROR;
} MPI_Status;

 
struct MPIR_T_enum_s;
struct MPIR_T_cvar_handle_s;
struct MPIR_T_pvar_handle_s;
struct MPIR_T_pvar_session_s;

typedef struct MPIR_T_enum_s * MPI_T_enum;
typedef struct MPIR_T_cvar_handle_s * MPI_T_cvar_handle;
typedef struct MPIR_T_pvar_handle_s * MPI_T_pvar_handle;
typedef struct MPIR_T_pvar_session_s * MPI_T_pvar_session;

 
extern  struct MPIR_T_pvar_handle_s * const MPI_T_PVAR_ALL_HANDLES ;



 
typedef enum MPIR_T_verbosity_t {
    
 
    MPIX_T_VERBOSITY_INVALID = 0,

     
    MPI_T_VERBOSITY_USER_BASIC = 221,
    MPI_T_VERBOSITY_USER_DETAIL,
    MPI_T_VERBOSITY_USER_ALL,

    MPI_T_VERBOSITY_TUNER_BASIC,
    MPI_T_VERBOSITY_TUNER_DETAIL,
    MPI_T_VERBOSITY_TUNER_ALL,

    MPI_T_VERBOSITY_MPIDEV_BASIC,
    MPI_T_VERBOSITY_MPIDEV_DETAIL,
    MPI_T_VERBOSITY_MPIDEV_ALL
} MPIR_T_verbosity_t;

typedef enum MPIR_T_bind_t {
    
 
    MPIX_T_BIND_INVALID = 0,

     
    MPI_T_BIND_NO_OBJECT = 9700,
    MPI_T_BIND_MPI_COMM,
    MPI_T_BIND_MPI_DATATYPE,
    MPI_T_BIND_MPI_ERRHANDLER,
    MPI_T_BIND_MPI_FILE,
    MPI_T_BIND_MPI_GROUP,
    MPI_T_BIND_MPI_OP,
    MPI_T_BIND_MPI_REQUEST,
    MPI_T_BIND_MPI_WIN,
    MPI_T_BIND_MPI_MESSAGE,
    MPI_T_BIND_MPI_INFO
} MPIR_T_bind_t;

typedef enum MPIR_T_scope_t {
    
 
    MPIX_T_SCOPE_INVALID = 0,

     
    MPI_T_SCOPE_CONSTANT = 60438,
    MPI_T_SCOPE_READONLY,
    MPI_T_SCOPE_LOCAL,
    MPI_T_SCOPE_GROUP,
    MPI_T_SCOPE_GROUP_EQ,
    MPI_T_SCOPE_ALL,
    MPI_T_SCOPE_ALL_EQ
} MPIR_T_scope_t;

typedef enum MPIR_T_pvar_class_t {
    
 
    MPIX_T_PVAR_CLASS_INVALID = 0,

     
    MPIR_T_PVAR_CLASS_FIRST = 240,
    MPI_T_PVAR_CLASS_STATE = MPIR_T_PVAR_CLASS_FIRST,
    MPI_T_PVAR_CLASS_LEVEL,
    MPI_T_PVAR_CLASS_SIZE,
    MPI_T_PVAR_CLASS_PERCENTAGE,
    MPI_T_PVAR_CLASS_HIGHWATERMARK,
    MPI_T_PVAR_CLASS_LOWWATERMARK,
    MPI_T_PVAR_CLASS_COUNTER,
    MPI_T_PVAR_CLASS_AGGREGATE,
    MPI_T_PVAR_CLASS_TIMER,
    MPI_T_PVAR_CLASS_GENERIC,
    MPIR_T_PVAR_CLASS_LAST,
    MPIR_T_PVAR_CLASS_NUMBER = MPIR_T_PVAR_CLASS_LAST - MPIR_T_PVAR_CLASS_FIRST
} MPIR_T_pvar_class_t;

 

 

 


 
extern  MPI_Fint * MPI_F_STATUS_IGNORE ;
extern  MPI_Fint * MPI_F_STATUSES_IGNORE ;



 


 



 
typedef struct {
    MPI_Fint count_lo;
    MPI_Fint count_hi_and_cancelled;
    MPI_Fint MPI_SOURCE;
    MPI_Fint MPI_TAG;
    MPI_Fint MPI_ERROR;
} MPI_F08_status;

extern  MPI_F08_status MPIR_F08_MPI_STATUS_IGNORE_OBJ ;
extern  MPI_F08_status MPIR_F08_MPI_STATUSES_IGNORE_OBJ[1] ;
extern  int MPIR_F08_MPI_IN_PLACE ;
extern  int MPIR_F08_MPI_BOTTOM ;

 
extern  MPI_F08_status *MPI_F08_STATUS_IGNORE ;
extern  MPI_F08_status *MPI_F08_STATUSES_IGNORE ;

 

 
typedef int (MPI_Grequest_cancel_function)(void *, int); 
typedef int (MPI_Grequest_free_function)(void *); 
typedef int (MPI_Grequest_query_function)(void *, MPI_Status *); 
typedef int (MPIX_Grequest_poll_function)(void *, MPI_Status *);
typedef int (MPIX_Grequest_wait_function)(int, void **, double, MPI_Status *);

 
 

 

 

 

 


 

 


 




 







 

 
typedef int (MPI_Datarep_conversion_function)(void *, MPI_Datatype, int, 
             void *, MPI_Offset, void *);
typedef int (MPI_Datarep_extent_function)(MPI_Datatype datatype, MPI_Aint *,
                      void *);





 
 




 
 
 
int MPI_Send(const void *buf, int count, MPI_Datatype datatype, int dest, int tag,
             MPI_Comm comm)  ;
int MPI_Recv(void *buf, int count, MPI_Datatype datatype, int source, int tag,
             MPI_Comm comm, MPI_Status *status)  ;
int MPI_Get_count(const MPI_Status *status, MPI_Datatype datatype, int *count) ;
int MPI_Bsend(const void *buf, int count, MPI_Datatype datatype, int dest, int tag,
              MPI_Comm comm)  ;
int MPI_Ssend(const void *buf, int count, MPI_Datatype datatype, int dest, int tag,
              MPI_Comm comm)  ;
int MPI_Rsend(const void *buf, int count, MPI_Datatype datatype, int dest, int tag,
              MPI_Comm comm)  ;
int MPI_Buffer_attach(void *buffer, int size) ;
int MPI_Buffer_detach(void *buffer_addr, int *size) ;
int MPI_Isend(const void *buf, int count, MPI_Datatype datatype, int dest, int tag,
              MPI_Comm comm, MPI_Request *request)  ;
int MPI_Ibsend(const void *buf, int count, MPI_Datatype datatype, int dest, int tag,
               MPI_Comm comm, MPI_Request *request)  ;
int MPI_Issend(const void *buf, int count, MPI_Datatype datatype, int dest, int tag,
               MPI_Comm comm, MPI_Request *request)  ;
int MPI_Irsend(const void *buf, int count, MPI_Datatype datatype, int dest, int tag,
               MPI_Comm comm, MPI_Request *request)  ;
int MPI_Irecv(void *buf, int count, MPI_Datatype datatype, int source, int tag,
              MPI_Comm comm, MPI_Request *request)  ;
int MPI_Wait(MPI_Request *request, MPI_Status *status) ;
int MPI_Test(MPI_Request *request, int *flag, MPI_Status *status) ;
int MPI_Request_free(MPI_Request *request) ;
int MPI_Waitany(int count, MPI_Request array_of_requests[], int *indx, MPI_Status *status) ;
int MPI_Testany(int count, MPI_Request array_of_requests[], int *indx, int *flag,
                MPI_Status *status) ;
int MPI_Waitall(int count, MPI_Request array_of_requests[], MPI_Status array_of_statuses[]) ;
int MPI_Testall(int count, MPI_Request array_of_requests[], int *flag,
                MPI_Status array_of_statuses[]) ;
int MPI_Waitsome(int incount, MPI_Request array_of_requests[], int *outcount,
                 int array_of_indices[], MPI_Status array_of_statuses[]) ;
int MPI_Testsome(int incount, MPI_Request array_of_requests[], int *outcount,
                 int array_of_indices[], MPI_Status array_of_statuses[]) ;
int MPI_Iprobe(int source, int tag, MPI_Comm comm, int *flag, MPI_Status *status) ;
int MPI_Probe(int source, int tag, MPI_Comm comm, MPI_Status *status) ;
int MPI_Cancel(MPI_Request *request) ;
int MPI_Test_cancelled(const MPI_Status *status, int *flag) ;
int MPI_Send_init(const void *buf, int count, MPI_Datatype datatype, int dest, int tag,
                  MPI_Comm comm, MPI_Request *request)  ;
int MPI_Bsend_init(const void *buf, int count, MPI_Datatype datatype, int dest, int tag,
                   MPI_Comm comm, MPI_Request *request)  ;
int MPI_Ssend_init(const void *buf, int count, MPI_Datatype datatype, int dest, int tag,
                   MPI_Comm comm, MPI_Request *request)  ;
int MPI_Rsend_init(const void *buf, int count, MPI_Datatype datatype, int dest, int tag,
                   MPI_Comm comm, MPI_Request *request)  ;
int MPI_Recv_init(void *buf, int count, MPI_Datatype datatype, int source, int tag,
                  MPI_Comm comm, MPI_Request *request)  ;
int MPI_Start(MPI_Request *request) ;
int MPI_Startall(int count, MPI_Request array_of_requests[]) ;
int MPI_Sendrecv(const void *sendbuf, int sendcount, MPI_Datatype sendtype, int dest,
                 int sendtag, void *recvbuf, int recvcount, MPI_Datatype recvtype,
                 int source, int recvtag, MPI_Comm comm, MPI_Status *status)
                   ;
int MPI_Sendrecv_replace(void *buf, int count, MPI_Datatype datatype, int dest,
                         int sendtag, int source, int recvtag, MPI_Comm comm,
                         MPI_Status *status)  ;
int MPI_Type_contiguous(int count, MPI_Datatype oldtype, MPI_Datatype *newtype) ;
int MPI_Type_vector(int count, int blocklength, int stride, MPI_Datatype oldtype,
                    MPI_Datatype *newtype) ;
int MPI_Type_hvector(int count, int blocklength, MPI_Aint stride, MPI_Datatype oldtype,
                     MPI_Datatype *newtype) ;
int MPI_Type_indexed(int count, const int *array_of_blocklengths,
                     const int *array_of_displacements, MPI_Datatype oldtype,
                     MPI_Datatype *newtype) ;
int MPI_Type_hindexed(int count, int *array_of_blocklengths,
                      MPI_Aint *array_of_displacements, MPI_Datatype oldtype,
                      MPI_Datatype *newtype) ;
int MPI_Type_struct(int count, int *array_of_blocklengths,
                    MPI_Aint *array_of_displacements,
                    MPI_Datatype *array_of_types, MPI_Datatype *newtype) ;
int MPI_Address(void *location, MPI_Aint *address) ;
int MPI_Type_extent(MPI_Datatype datatype, MPI_Aint *extent) ;
int MPI_Type_size(MPI_Datatype datatype, int *size) ;
int MPI_Type_lb(MPI_Datatype datatype, MPI_Aint *displacement) ;
int MPI_Type_ub(MPI_Datatype datatype, MPI_Aint *displacement) ;
int MPI_Type_commit(MPI_Datatype *datatype) ;
int MPI_Type_free(MPI_Datatype *datatype) ;
int MPI_Get_elements(const MPI_Status *status, MPI_Datatype datatype, int *count) ;
int MPI_Pack(const void *inbuf, int incount, MPI_Datatype datatype, void *outbuf,
             int outsize, int *position, MPI_Comm comm)  ;
int MPI_Unpack(const void *inbuf, int insize, int *position, void *outbuf, int outcount,
               MPI_Datatype datatype, MPI_Comm comm)  ;
int MPI_Pack_size(int incount, MPI_Datatype datatype, MPI_Comm comm, int *size) ;
int MPI_Barrier(MPI_Comm comm) ;
int MPI_Bcast(void *buffer, int count, MPI_Datatype datatype, int root, MPI_Comm comm)
               ;
int MPI_Gather(const void *sendbuf, int sendcount, MPI_Datatype sendtype, void *recvbuf,
               int recvcount, MPI_Datatype recvtype, int root, MPI_Comm comm)
                 ;
int MPI_Gatherv(const void *sendbuf, int sendcount, MPI_Datatype sendtype, void *recvbuf,
                const int *recvcounts, const int *displs, MPI_Datatype recvtype, int root,
                MPI_Comm comm)
                  ;
int MPI_Scatter(const void *sendbuf, int sendcount, MPI_Datatype sendtype, void *recvbuf,
                int recvcount, MPI_Datatype recvtype, int root, MPI_Comm comm)
                  ;
int MPI_Scatterv(const void *sendbuf, const int *sendcounts, const int *displs,
                 MPI_Datatype sendtype, void *recvbuf, int recvcount, MPI_Datatype recvtype,
                 int root, MPI_Comm comm)
                   ;
int MPI_Allgather(const void *sendbuf, int sendcount, MPI_Datatype sendtype, void *recvbuf,
                  int recvcount, MPI_Datatype recvtype, MPI_Comm comm)
                    ;
int MPI_Allgatherv(const void *sendbuf, int sendcount, MPI_Datatype sendtype, void *recvbuf,
                   const int *recvcounts, const int *displs, MPI_Datatype recvtype, MPI_Comm comm)
                     ;
int MPI_Alltoall(const void *sendbuf, int sendcount, MPI_Datatype sendtype, void *recvbuf,
                 int recvcount, MPI_Datatype recvtype, MPI_Comm comm)
                   ;
int MPI_Alltoallv(const void *sendbuf, const int *sendcounts, const int *sdispls,
                  MPI_Datatype sendtype, void *recvbuf, const int *recvcounts,
                  const int *rdispls, MPI_Datatype recvtype, MPI_Comm comm)
                    ;
int MPI_Alltoallw(const void *sendbuf, const int sendcounts[], const int sdispls[],
                  const MPI_Datatype sendtypes[], void *recvbuf, const int recvcounts[],
                  const int rdispls[], const MPI_Datatype recvtypes[], MPI_Comm comm) ;
int MPI_Exscan(const void *sendbuf, void *recvbuf, int count, MPI_Datatype datatype,
               MPI_Op op, MPI_Comm comm)
                 ;
int MPI_Reduce(const void *sendbuf, void *recvbuf, int count, MPI_Datatype datatype,
               MPI_Op op, int root, MPI_Comm comm)
                 ;
int MPI_Op_create(MPI_User_function *user_fn, int commute, MPI_Op *op) ;
int MPI_Op_free(MPI_Op *op) ;
int MPI_Allreduce(const void *sendbuf, void *recvbuf, int count, MPI_Datatype datatype,
                  MPI_Op op, MPI_Comm comm)
                    ;
int MPI_Reduce_scatter(const void *sendbuf, void *recvbuf, const int recvcounts[],
                       MPI_Datatype datatype, MPI_Op op, MPI_Comm comm)
                         ;
int MPI_Scan(const void *sendbuf, void *recvbuf, int count, MPI_Datatype datatype, MPI_Op op,
             MPI_Comm comm)
               ;
int MPI_Group_size(MPI_Group group, int *size) ;
int MPI_Group_rank(MPI_Group group, int *rank) ;
int MPI_Group_translate_ranks(MPI_Group group1, int n, const int ranks1[], MPI_Group group2,
                              int ranks2[]) ;
int MPI_Group_compare(MPI_Group group1, MPI_Group group2, int *result) ;
int MPI_Comm_group(MPI_Comm comm, MPI_Group *group) ;
int MPI_Group_union(MPI_Group group1, MPI_Group group2, MPI_Group *newgroup) ;
int MPI_Group_intersection(MPI_Group group1, MPI_Group group2, MPI_Group *newgroup) ;
int MPI_Group_difference(MPI_Group group1, MPI_Group group2, MPI_Group *newgroup) ;
int MPI_Group_incl(MPI_Group group, int n, const int ranks[], MPI_Group *newgroup) ;
int MPI_Group_excl(MPI_Group group, int n, const int ranks[], MPI_Group *newgroup) ;
int MPI_Group_range_incl(MPI_Group group, int n, int ranges[][3], MPI_Group *newgroup) ;
int MPI_Group_range_excl(MPI_Group group, int n, int ranges[][3], MPI_Group *newgroup) ;
int MPI_Group_free(MPI_Group *group) ;
int MPI_Comm_size(MPI_Comm comm, int *size) ;
int MPI_Comm_rank(MPI_Comm comm, int *rank) ;
int MPI_Comm_compare(MPI_Comm comm1, MPI_Comm comm2, int *result) ;
int MPI_Comm_dup(MPI_Comm comm, MPI_Comm *newcomm) ;
int MPI_Comm_dup_with_info(MPI_Comm comm, MPI_Info info, MPI_Comm *newcomm) ;
int MPI_Comm_create(MPI_Comm comm, MPI_Group group, MPI_Comm *newcomm) ;
int MPI_Comm_split(MPI_Comm comm, int color, int key, MPI_Comm *newcomm) ;
int MPI_Comm_free(MPI_Comm *comm) ;
int MPI_Comm_test_inter(MPI_Comm comm, int *flag) ;
int MPI_Comm_remote_size(MPI_Comm comm, int *size) ;
int MPI_Comm_remote_group(MPI_Comm comm, MPI_Group *group) ;
int MPI_Intercomm_create(MPI_Comm local_comm, int local_leader, MPI_Comm peer_comm,
                         int remote_leader, int tag, MPI_Comm *newintercomm) ;
int MPI_Intercomm_merge(MPI_Comm intercomm, int high, MPI_Comm *newintracomm) ;
int MPI_Keyval_create(MPI_Copy_function *copy_fn, MPI_Delete_function *delete_fn,
                      int *keyval, void *extra_state) ;
int MPI_Keyval_free(int *keyval) ;
int MPI_Attr_put(MPI_Comm comm, int keyval, void *attribute_val) ;
int MPI_Attr_get(MPI_Comm comm, int keyval, void *attribute_val, int *flag) ;
int MPI_Attr_delete(MPI_Comm comm, int keyval) ;
int MPI_Topo_test(MPI_Comm comm, int *status) ;
int MPI_Cart_create(MPI_Comm comm_old, int ndims, const int dims[], const int periods[],
                    int reorder, MPI_Comm *comm_cart) ;
int MPI_Dims_create(int nnodes, int ndims, int dims[]) ;
int MPI_Graph_create(MPI_Comm comm_old, int nnodes, const int indx[], const int edges[],
                     int reorder, MPI_Comm *comm_graph) ;
int MPI_Graphdims_get(MPI_Comm comm, int *nnodes, int *nedges) ;
int MPI_Graph_get(MPI_Comm comm, int maxindex, int maxedges, int indx[], int edges[]) ;
int MPI_Cartdim_get(MPI_Comm comm, int *ndims) ;
int MPI_Cart_get(MPI_Comm comm, int maxdims, int dims[], int periods[], int coords[]) ;
int MPI_Cart_rank(MPI_Comm comm, const int coords[], int *rank) ;
int MPI_Cart_coords(MPI_Comm comm, int rank, int maxdims, int coords[]) ;
int MPI_Graph_neighbors_count(MPI_Comm comm, int rank, int *nneighbors) ;
int MPI_Graph_neighbors(MPI_Comm comm, int rank, int maxneighbors, int neighbors[]) ;
int MPI_Cart_shift(MPI_Comm comm, int direction, int disp, int *rank_source, int *rank_dest) ;
int MPI_Cart_sub(MPI_Comm comm, const int remain_dims[], MPI_Comm *newcomm) ;
int MPI_Cart_map(MPI_Comm comm, int ndims, const int dims[], const int periods[], int *newrank) ;
int MPI_Graph_map(MPI_Comm comm, int nnodes, const int indx[], const int edges[], int *newrank) ;
int MPI_Get_processor_name(char *name, int *resultlen) ;
int MPI_Get_version(int *version, int *subversion) ;
int MPI_Get_library_version(char *version, int *resultlen) ;
int MPI_Errhandler_create(MPI_Handler_function *function, MPI_Errhandler *errhandler) ;
int MPI_Errhandler_set(MPI_Comm comm, MPI_Errhandler errhandler) ;
int MPI_Errhandler_get(MPI_Comm comm, MPI_Errhandler *errhandler) ;
int MPI_Errhandler_free(MPI_Errhandler *errhandler) ;
int MPI_Error_string(int errorcode, char *string, int *resultlen) ;
int MPI_Error_class(int errorcode, int *errorclass) ;
double MPI_Wtime(void) ;
double MPI_Wtick(void) ;
int MPI_Init(int *argc, char ***argv) ;
int MPI_Finalize(void) ;
int MPI_Initialized(int *flag) ;
int MPI_Abort(MPI_Comm comm, int errorcode) ;


 
int MPI_Pcontrol(const int level, ...) ;
int MPIR_Dup_fn(MPI_Comm oldcomm, int keyval, void *extra_state, void *attribute_val_in,
               void *attribute_val_out, int *flag) ;

 
int MPI_Close_port(const char *port_name) ;
int MPI_Comm_accept(const char *port_name, MPI_Info info, int root, MPI_Comm comm,
                    MPI_Comm *newcomm) ;
int MPI_Comm_connect(const char *port_name, MPI_Info info, int root, MPI_Comm comm,
                     MPI_Comm *newcomm) ;
int MPI_Comm_disconnect(MPI_Comm *comm) ;
int MPI_Comm_get_parent(MPI_Comm *parent) ;
int MPI_Comm_join(int fd, MPI_Comm *intercomm) ;
int MPI_Comm_spawn(const char *command, char *argv[], int maxprocs, MPI_Info info, int root,
                   MPI_Comm comm, MPI_Comm *intercomm, int array_of_errcodes[]) ;
int MPI_Comm_spawn_multiple(int count, char *array_of_commands[], char **array_of_argv[],
                            const int array_of_maxprocs[], const MPI_Info array_of_info[],
                            int root, MPI_Comm comm, MPI_Comm *intercomm, int array_of_errcodes[]) ;
int MPI_Lookup_name(const char *service_name, MPI_Info info, char *port_name) ;
int MPI_Open_port(MPI_Info info, char *port_name) ;
int MPI_Publish_name(const char *service_name, MPI_Info info, const char *port_name) ;
int MPI_Unpublish_name(const char *service_name, MPI_Info info, const char *port_name) ;
int MPI_Comm_set_info(MPI_Comm comm, MPI_Info info) ;
int MPI_Comm_get_info(MPI_Comm comm, MPI_Info *info) ;

 
int MPI_Accumulate(const void *origin_addr, int origin_count, MPI_Datatype origin_datatype,
                   int target_rank, MPI_Aint target_disp, int target_count,
                   MPI_Datatype target_datatype, MPI_Op op, MPI_Win win)
                    ;
int MPI_Get(void *origin_addr, int origin_count, MPI_Datatype origin_datatype,
            int target_rank, MPI_Aint target_disp, int target_count,
            MPI_Datatype target_datatype, MPI_Win win)  ;
int MPI_Put(const void *origin_addr, int origin_count, MPI_Datatype origin_datatype,
            int target_rank, MPI_Aint target_disp, int target_count,
            MPI_Datatype target_datatype, MPI_Win win)  ;
int MPI_Win_complete(MPI_Win win) ;
int MPI_Win_create(void *base, MPI_Aint size, int disp_unit, MPI_Info info, MPI_Comm comm,
                   MPI_Win *win) ;
int MPI_Win_fence(int assert, MPI_Win win) ;
int MPI_Win_free(MPI_Win *win) ;
int MPI_Win_get_group(MPI_Win win, MPI_Group *group) ;
int MPI_Win_lock(int lock_type, int rank, int assert, MPI_Win win) ;
int MPI_Win_post(MPI_Group group, int assert, MPI_Win win) ;
int MPI_Win_start(MPI_Group group, int assert, MPI_Win win) ;
int MPI_Win_test(MPI_Win win, int *flag) ;
int MPI_Win_unlock(int rank, MPI_Win win) ;
int MPI_Win_wait(MPI_Win win) ;

 
int MPI_Win_allocate(MPI_Aint size, int disp_unit, MPI_Info info, MPI_Comm comm, void *baseptr,
                     MPI_Win *win) ;
int MPI_Win_allocate_shared(MPI_Aint size, int disp_unit, MPI_Info info, MPI_Comm comm,
                            void *baseptr, MPI_Win *win) ;
int MPI_Win_shared_query(MPI_Win win, int rank, MPI_Aint *size, int *disp_unit, void *baseptr) ;
int MPI_Win_create_dynamic(MPI_Info info, MPI_Comm comm, MPI_Win *win) ;
int MPI_Win_attach(MPI_Win win, void *base, MPI_Aint size) ;
int MPI_Win_detach(MPI_Win win, const void *base) ;
int MPI_Win_get_info(MPI_Win win, MPI_Info *info_used) ;
int MPI_Win_set_info(MPI_Win win, MPI_Info info) ;
int MPI_Get_accumulate(const void *origin_addr, int origin_count,
                        MPI_Datatype origin_datatype, void *result_addr, int result_count,
                        MPI_Datatype result_datatype, int target_rank, MPI_Aint target_disp,
                        int target_count, MPI_Datatype target_datatype, MPI_Op op, MPI_Win win)
                        
                         ;
int MPI_Fetch_and_op(const void *origin_addr, void *result_addr,
                      MPI_Datatype datatype, int target_rank, MPI_Aint target_disp,
                      MPI_Op op, MPI_Win win)
                       ;
int MPI_Compare_and_swap(const void *origin_addr, const void *compare_addr,
                          void *result_addr, MPI_Datatype datatype, int target_rank,
                          MPI_Aint target_disp, MPI_Win win)
                          
                          
                           ;
int MPI_Rput(const void *origin_addr, int origin_count,
              MPI_Datatype origin_datatype, int target_rank, MPI_Aint target_disp,
              int target_count, MPI_Datatype target_datatype, MPI_Win win,
              MPI_Request *request)
               ;
int MPI_Rget(void *origin_addr, int origin_count,
              MPI_Datatype origin_datatype, int target_rank, MPI_Aint target_disp,
              int target_count, MPI_Datatype target_datatype, MPI_Win win,
              MPI_Request *request)
               ;
int MPI_Raccumulate(const void *origin_addr, int origin_count,
                     MPI_Datatype origin_datatype, int target_rank, MPI_Aint target_disp,
                     int target_count, MPI_Datatype target_datatype, MPI_Op op, MPI_Win win,
                     MPI_Request *request)
                      ;
int MPI_Rget_accumulate(const void *origin_addr, int origin_count,
                         MPI_Datatype origin_datatype, void *result_addr, int result_count,
                         MPI_Datatype result_datatype, int target_rank, MPI_Aint target_disp,
                         int target_count, MPI_Datatype target_datatype, MPI_Op op, MPI_Win win,
                         MPI_Request *request)
                         
                          ;
int MPI_Win_lock_all(int assert, MPI_Win win) ;
int MPI_Win_unlock_all(MPI_Win win) ;
int MPI_Win_flush(int rank, MPI_Win win) ;
int MPI_Win_flush_all(MPI_Win win) ;
int MPI_Win_flush_local(int rank, MPI_Win win) ;
int MPI_Win_flush_local_all(MPI_Win win) ;
int MPI_Win_sync(MPI_Win win) ;
 
 
int MPI_Add_error_class(int *errorclass) ;
int MPI_Add_error_code(int errorclass, int *errorcode) ;
int MPI_Add_error_string(int errorcode, const char *string) ;
int MPI_Comm_call_errhandler(MPI_Comm comm, int errorcode) ;
int MPI_Comm_create_keyval(MPI_Comm_copy_attr_function *comm_copy_attr_fn,
                           MPI_Comm_delete_attr_function *comm_delete_attr_fn, int *comm_keyval,
                           void *extra_state) ;
int MPI_Comm_delete_attr(MPI_Comm comm, int comm_keyval) ;
int MPI_Comm_free_keyval(int *comm_keyval) ;
int MPI_Comm_get_attr(MPI_Comm comm, int comm_keyval, void *attribute_val, int *flag) ;
int MPI_Comm_get_name(MPI_Comm comm, char *comm_name, int *resultlen) ;
int MPI_Comm_set_attr(MPI_Comm comm, int comm_keyval, void *attribute_val) ;
int MPI_Comm_set_name(MPI_Comm comm, const char *comm_name) ;
int MPI_File_call_errhandler(MPI_File fh, int errorcode) ;
int MPI_Grequest_complete(MPI_Request request) ;
int MPI_Grequest_start(MPI_Grequest_query_function *query_fn, MPI_Grequest_free_function *free_fn,
                       MPI_Grequest_cancel_function *cancel_fn, void *extra_state,
                       MPI_Request *request) ;
int MPI_Init_thread(int *argc, char ***argv, int required, int *provided) ;
int MPI_Is_thread_main(int *flag) ;
int MPI_Query_thread(int *provided) ;
int MPI_Status_set_cancelled(MPI_Status *status, int flag) ;
int MPI_Status_set_elements(MPI_Status *status, MPI_Datatype datatype, int count) ;
int MPI_Type_create_keyval(MPI_Type_copy_attr_function *type_copy_attr_fn,
                           MPI_Type_delete_attr_function *type_delete_attr_fn,
                           int *type_keyval, void *extra_state) ;
int MPI_Type_delete_attr(MPI_Datatype datatype, int type_keyval) ;
int MPI_Type_dup(MPI_Datatype oldtype, MPI_Datatype *newtype) ;
int MPI_Type_free_keyval(int *type_keyval) ;
int MPI_Type_get_attr(MPI_Datatype datatype, int type_keyval, void *attribute_val, int *flag) ;
int MPI_Type_get_contents(MPI_Datatype datatype, int max_integers, int max_addresses,
                          int max_datatypes, int array_of_integers[],
                          MPI_Aint array_of_addresses[], MPI_Datatype array_of_datatypes[]) ;
int MPI_Type_get_envelope(MPI_Datatype datatype, int *num_integers, int *num_addresses,
                          int *num_datatypes, int *combiner) ;
int MPI_Type_get_name(MPI_Datatype datatype, char *type_name, int *resultlen) ;
int MPI_Type_set_attr(MPI_Datatype datatype, int type_keyval, void *attribute_val) ;
int MPI_Type_set_name(MPI_Datatype datatype, const char *type_name) ;
int MPI_Type_match_size(int typeclass, int size, MPI_Datatype *datatype) ;
int MPI_Win_call_errhandler(MPI_Win win, int errorcode) ;
int MPI_Win_create_keyval(MPI_Win_copy_attr_function *win_copy_attr_fn,
                          MPI_Win_delete_attr_function *win_delete_attr_fn, int *win_keyval,
                          void *extra_state) ;
int MPI_Win_delete_attr(MPI_Win win, int win_keyval) ;
int MPI_Win_free_keyval(int *win_keyval) ;
int MPI_Win_get_attr(MPI_Win win, int win_keyval, void *attribute_val, int *flag) ;
int MPI_Win_get_name(MPI_Win win, char *win_name, int *resultlen) ;
int MPI_Win_set_attr(MPI_Win win, int win_keyval, void *attribute_val) ;
int MPI_Win_set_name(MPI_Win win, const char *win_name) ;

int MPI_Alloc_mem(MPI_Aint size, MPI_Info info, void *baseptr) ;
int MPI_Comm_create_errhandler(MPI_Comm_errhandler_function *comm_errhandler_fn,
                               MPI_Errhandler *errhandler) ;
int MPI_Comm_get_errhandler(MPI_Comm comm, MPI_Errhandler *errhandler) ;
int MPI_Comm_set_errhandler(MPI_Comm comm, MPI_Errhandler errhandler) ;
int MPI_File_create_errhandler(MPI_File_errhandler_function *file_errhandler_fn,
                               MPI_Errhandler *errhandler) ;
int MPI_File_get_errhandler(MPI_File file, MPI_Errhandler *errhandler) ;
int MPI_File_set_errhandler(MPI_File file, MPI_Errhandler errhandler) ;
int MPI_Finalized(int *flag) ;
int MPI_Free_mem(void *base) ;
int MPI_Get_address(const void *location, MPI_Aint *address) ;
int MPI_Info_create(MPI_Info *info) ;
int MPI_Info_delete(MPI_Info info, const char *key) ;
int MPI_Info_dup(MPI_Info info, MPI_Info *newinfo) ;
int MPI_Info_free(MPI_Info *info) ;
int MPI_Info_get(MPI_Info info, const char *key, int valuelen, char *value, int *flag) ;
int MPI_Info_get_nkeys(MPI_Info info, int *nkeys) ;
int MPI_Info_get_nthkey(MPI_Info info, int n, char *key) ;
int MPI_Info_get_valuelen(MPI_Info info, const char *key, int *valuelen, int *flag) ;
int MPI_Info_set(MPI_Info info, const char *key, const char *value) ;
int MPI_Pack_external(const char datarep[], const void *inbuf, int incount,
                      MPI_Datatype datatype, void *outbuf, MPI_Aint outsize, MPI_Aint *position)
                       ;
int MPI_Pack_external_size(const char datarep[], int incount, MPI_Datatype datatype,
                           MPI_Aint *size) ;
int MPI_Request_get_status(MPI_Request request, int *flag, MPI_Status *status) ;
int MPI_Status_c2f(const MPI_Status *c_status, MPI_Fint *f_status) ;
int MPI_Status_f2c(const MPI_Fint *f_status, MPI_Status *c_status) ;
int MPI_Type_create_darray(int size, int rank, int ndims, const int array_of_gsizes[],
                           const int array_of_distribs[], const int array_of_dargs[],
                           const int array_of_psizes[], int order, MPI_Datatype oldtype,
                           MPI_Datatype *newtype) ;
int MPI_Type_create_hindexed(int count, const int array_of_blocklengths[],
                             const MPI_Aint array_of_displacements[], MPI_Datatype oldtype,
                             MPI_Datatype *newtype) ;
int MPI_Type_create_hvector(int count, int blocklength, MPI_Aint stride, MPI_Datatype oldtype,
                            MPI_Datatype *newtype) ;
int MPI_Type_create_indexed_block(int count, int blocklength, const int array_of_displacements[],
                                  MPI_Datatype oldtype, MPI_Datatype *newtype) ;
int MPI_Type_create_hindexed_block(int count, int blocklength,
                                   const MPI_Aint array_of_displacements[],
                                   MPI_Datatype oldtype, MPI_Datatype *newtype) ;
int MPI_Type_create_resized(MPI_Datatype oldtype, MPI_Aint lb, MPI_Aint extent,
                            MPI_Datatype *newtype) ;
int MPI_Type_create_struct(int count, const int array_of_blocklengths[],
                           const MPI_Aint array_of_displacements[],
                           const MPI_Datatype array_of_types[], MPI_Datatype *newtype) ;
int MPI_Type_create_subarray(int ndims, const int array_of_sizes[],
                             const int array_of_subsizes[], const int array_of_starts[],
                             int order, MPI_Datatype oldtype, MPI_Datatype *newtype) ;
int MPI_Type_get_extent(MPI_Datatype datatype, MPI_Aint *lb, MPI_Aint *extent) ;
int MPI_Type_get_true_extent(MPI_Datatype datatype, MPI_Aint *true_lb, MPI_Aint *true_extent) ;
int MPI_Unpack_external(const char datarep[], const void *inbuf, MPI_Aint insize,
                        MPI_Aint *position, void *outbuf, int outcount, MPI_Datatype datatype)
                         ;
int MPI_Win_create_errhandler(MPI_Win_errhandler_function *win_errhandler_fn,
                              MPI_Errhandler *errhandler) ;
int MPI_Win_get_errhandler(MPI_Win win, MPI_Errhandler *errhandler) ;
int MPI_Win_set_errhandler(MPI_Win win, MPI_Errhandler errhandler) ;



 
int MPI_Type_create_f90_integer(int range, MPI_Datatype *newtype) ;
int MPI_Type_create_f90_real(int precision, int range, MPI_Datatype *newtype) ;
int MPI_Type_create_f90_complex(int precision, int range, MPI_Datatype *newtype) ;

int MPI_Reduce_local(const void *inbuf, void *inoutbuf, int count, MPI_Datatype datatype,
                     MPI_Op op)
                       ;
int MPI_Op_commutative(MPI_Op op, int *commute) ;
int MPI_Reduce_scatter_block(const void *sendbuf, void *recvbuf, int recvcount,
                             MPI_Datatype datatype, MPI_Op op, MPI_Comm comm)
                             
                              ;
int MPI_Dist_graph_create_adjacent(MPI_Comm comm_old, int indegree, const int sources[],
                                   const int sourceweights[], int outdegree,
                                   const int destinations[], const int destweights[],
                                   MPI_Info info, int reorder, MPI_Comm *comm_dist_graph) ;
int MPI_Dist_graph_create(MPI_Comm comm_old, int n, const int sources[], const int degrees[],
                          const int destinations[], const int weights[], MPI_Info info,
                          int reorder, MPI_Comm *comm_dist_graph) ;
int MPI_Dist_graph_neighbors_count(MPI_Comm comm, int *indegree, int *outdegree, int *weighted) ;
int MPI_Dist_graph_neighbors(MPI_Comm comm, int maxindegree, int sources[], int sourceweights[],
                             int maxoutdegree, int destinations[], int destweights[]) ;

 
int MPI_Improbe(int source, int tag, MPI_Comm comm, int *flag, MPI_Message *message,
                MPI_Status *status) ;
int MPI_Imrecv(void *buf, int count, MPI_Datatype datatype, MPI_Message *message,
               MPI_Request *request)  ;
int MPI_Mprobe(int source, int tag, MPI_Comm comm, MPI_Message *message, MPI_Status *status) ;
int MPI_Mrecv(void *buf, int count, MPI_Datatype datatype, MPI_Message *message,
              MPI_Status *status)  ;

 
int MPI_Comm_idup(MPI_Comm comm, MPI_Comm *newcomm, MPI_Request *request) ;
int MPI_Ibarrier(MPI_Comm comm, MPI_Request *request) ;
int MPI_Ibcast(void *buffer, int count, MPI_Datatype datatype, int root, MPI_Comm comm,
               MPI_Request *request)  ;
int MPI_Igather(const void *sendbuf, int sendcount, MPI_Datatype sendtype, void *recvbuf,
                int recvcount, MPI_Datatype recvtype, int root, MPI_Comm comm,
                MPI_Request *request)
                  ;
int MPI_Igatherv(const void *sendbuf, int sendcount, MPI_Datatype sendtype, void *recvbuf,
                 const int recvcounts[], const int displs[], MPI_Datatype recvtype, int root,
                 MPI_Comm comm, MPI_Request *request)
                   ;
int MPI_Iscatter(const void *sendbuf, int sendcount, MPI_Datatype sendtype, void *recvbuf,
                 int recvcount, MPI_Datatype recvtype, int root, MPI_Comm comm,
                 MPI_Request *request)
                   ;
int MPI_Iscatterv(const void *sendbuf, const int sendcounts[], const int displs[],
                  MPI_Datatype sendtype, void *recvbuf, int recvcount, MPI_Datatype recvtype,
                  int root, MPI_Comm comm, MPI_Request *request)
                    ;
int MPI_Iallgather(const void *sendbuf, int sendcount, MPI_Datatype sendtype, void *recvbuf,
                   int recvcount, MPI_Datatype recvtype, MPI_Comm comm, MPI_Request *request)
                     ;
int MPI_Iallgatherv(const void *sendbuf, int sendcount, MPI_Datatype sendtype, void *recvbuf,
                    const int recvcounts[], const int displs[], MPI_Datatype recvtype,
                    MPI_Comm comm, MPI_Request *request)
                      ;
int MPI_Ialltoall(const void *sendbuf, int sendcount, MPI_Datatype sendtype, void *recvbuf,
                  int recvcount, MPI_Datatype recvtype, MPI_Comm comm, MPI_Request *request)
                    ;
int MPI_Ialltoallv(const void *sendbuf, const int sendcounts[], const int sdispls[],
                   MPI_Datatype sendtype, void *recvbuf, const int recvcounts[],
                   const int rdispls[], MPI_Datatype recvtype, MPI_Comm comm,
                   MPI_Request *request)
                     ;
int MPI_Ialltoallw(const void *sendbuf, const int sendcounts[], const int sdispls[],
                   const MPI_Datatype sendtypes[], void *recvbuf, const int recvcounts[],
                   const int rdispls[], const MPI_Datatype recvtypes[], MPI_Comm comm,
                   MPI_Request *request) ;
int MPI_Ireduce(const void *sendbuf, void *recvbuf, int count, MPI_Datatype datatype,
                MPI_Op op, int root, MPI_Comm comm, MPI_Request *request)
                  ;
int MPI_Iallreduce(const void *sendbuf, void *recvbuf, int count, MPI_Datatype datatype,
                   MPI_Op op, MPI_Comm comm, MPI_Request *request)
                     ;
int MPI_Ireduce_scatter(const void *sendbuf, void *recvbuf, const int recvcounts[],
                        MPI_Datatype datatype, MPI_Op op, MPI_Comm comm, MPI_Request *request)
                          ;
int MPI_Ireduce_scatter_block(const void *sendbuf, void *recvbuf, int recvcount,
                              MPI_Datatype datatype, MPI_Op op, MPI_Comm comm,
                              MPI_Request *request)
                              
                               ;
int MPI_Iscan(const void *sendbuf, void *recvbuf, int count, MPI_Datatype datatype, MPI_Op op,
              MPI_Comm comm, MPI_Request *request)
                ;
int MPI_Iexscan(const void *sendbuf, void *recvbuf, int count, MPI_Datatype datatype,
                MPI_Op op, MPI_Comm comm, MPI_Request *request)
                  ;

 
int MPI_Ineighbor_allgather(const void *sendbuf, int sendcount, MPI_Datatype sendtype,
                            void *recvbuf, int recvcount, MPI_Datatype recvtype,
                            MPI_Comm comm, MPI_Request *request)
                            
                             ;
int MPI_Ineighbor_allgatherv(const void *sendbuf, int sendcount, MPI_Datatype sendtype,
                             void *recvbuf, const int recvcounts[], const int displs[],
                             MPI_Datatype recvtype, MPI_Comm comm, MPI_Request *request)
                             
                              ;
int MPI_Ineighbor_alltoall(const void *sendbuf, int sendcount, MPI_Datatype sendtype,
                           void *recvbuf, int recvcount, MPI_Datatype recvtype, MPI_Comm comm,
                           MPI_Request *request)
                           
                            ;
int MPI_Ineighbor_alltoallv(const void *sendbuf, const int sendcounts[], const int sdispls[],
                            MPI_Datatype sendtype, void *recvbuf, const int recvcounts[],
                            const int rdispls[], MPI_Datatype recvtype, MPI_Comm comm,
                            MPI_Request *request)
                            
                             ;
int MPI_Ineighbor_alltoallw(const void *sendbuf, const int sendcounts[],
                            const MPI_Aint sdispls[], const MPI_Datatype sendtypes[],
                            void *recvbuf, const int recvcounts[], const MPI_Aint rdispls[],
                            const MPI_Datatype recvtypes[], MPI_Comm comm, MPI_Request *request) ;
int MPI_Neighbor_allgather(const void *sendbuf, int sendcount, MPI_Datatype sendtype,
                           void *recvbuf, int recvcount, MPI_Datatype recvtype, MPI_Comm comm)
                           
                            ;
int MPI_Neighbor_allgatherv(const void *sendbuf, int sendcount, MPI_Datatype sendtype,
                            void *recvbuf, const int recvcounts[], const int displs[],
                            MPI_Datatype recvtype, MPI_Comm comm)
                            
                             ;
int MPI_Neighbor_alltoall(const void *sendbuf, int sendcount, MPI_Datatype sendtype,
                          void *recvbuf, int recvcount, MPI_Datatype recvtype, MPI_Comm comm)
                          
                           ;
int MPI_Neighbor_alltoallv(const void *sendbuf, const int sendcounts[], const int sdispls[],
                           MPI_Datatype sendtype, void *recvbuf, const int recvcounts[],
                           const int rdispls[], MPI_Datatype recvtype, MPI_Comm comm)
                           
                            ;
int MPI_Neighbor_alltoallw(const void *sendbuf, const int sendcounts[], const MPI_Aint sdispls[],
                           const MPI_Datatype sendtypes[], void *recvbuf, const int recvcounts[],
                           const MPI_Aint rdispls[], const MPI_Datatype recvtypes[], MPI_Comm comm) ;

 
int MPI_Comm_split_type(MPI_Comm comm, int split_type, int key, MPI_Info info, MPI_Comm *newcomm) ;

 
int MPI_Get_elements_x(const MPI_Status *status, MPI_Datatype datatype, MPI_Count *count) ;
int MPI_Status_set_elements_x(MPI_Status *status, MPI_Datatype datatype, MPI_Count count) ;
int MPI_Type_get_extent_x(MPI_Datatype datatype, MPI_Count *lb, MPI_Count *extent) ;
int MPI_Type_get_true_extent_x(MPI_Datatype datatype, MPI_Count *lb, MPI_Count *extent) ;
int MPI_Type_size_x(MPI_Datatype datatype, MPI_Count *size) ;

 
int MPI_Comm_create_group(MPI_Comm comm, MPI_Group group, int tag, MPI_Comm *newcomm) ;

 
MPI_Aint MPI_Aint_add(MPI_Aint base, MPI_Aint disp) ;
MPI_Aint MPI_Aint_diff(MPI_Aint addr1, MPI_Aint addr2) ;

 

 
 
int MPI_T_init_thread(int required, int *provided) ;
int MPI_T_finalize(void) ;
int MPI_T_enum_get_info(MPI_T_enum enumtype, int *num, char *name, int *name_len) ;
int MPI_T_enum_get_item(MPI_T_enum enumtype, int indx, int *value, char *name, int *name_len) ;
int MPI_T_cvar_get_num(int *num_cvar) ;
int MPI_T_cvar_get_info(int cvar_index, char *name, int *name_len, int *verbosity,
                        MPI_Datatype *datatype, MPI_T_enum *enumtype, char *desc, int *desc_len,
                        int *binding, int *scope) ;
int MPI_T_cvar_handle_alloc(int cvar_index, void *obj_handle, MPI_T_cvar_handle *handle,
                            int *count) ;
int MPI_T_cvar_handle_free(MPI_T_cvar_handle *handle) ;
int MPI_T_cvar_read(MPI_T_cvar_handle handle, void *buf) ;
int MPI_T_cvar_write(MPI_T_cvar_handle handle, const void *buf) ;
int MPI_T_pvar_get_num(int *num_pvar) ;
int MPI_T_pvar_get_info(int pvar_index, char *name, int *name_len, int *verbosity, int *var_class,
                        MPI_Datatype *datatype, MPI_T_enum *enumtype, char *desc, int *desc_len,
                        int *binding, int *readonly, int *continuous, int *atomic) ;
int MPI_T_pvar_session_create(MPI_T_pvar_session *session) ;
int MPI_T_pvar_session_free(MPI_T_pvar_session *session) ;
int MPI_T_pvar_handle_alloc(MPI_T_pvar_session session, int pvar_index, void *obj_handle,
                            MPI_T_pvar_handle *handle, int *count) ;
int MPI_T_pvar_handle_free(MPI_T_pvar_session session, MPI_T_pvar_handle *handle) ;
int MPI_T_pvar_start(MPI_T_pvar_session session, MPI_T_pvar_handle handle) ;
int MPI_T_pvar_stop(MPI_T_pvar_session session, MPI_T_pvar_handle handle) ;
int MPI_T_pvar_read(MPI_T_pvar_session session, MPI_T_pvar_handle handle, void *buf) ;
int MPI_T_pvar_write(MPI_T_pvar_session session, MPI_T_pvar_handle handle, const void *buf) ;
int MPI_T_pvar_reset(MPI_T_pvar_session session, MPI_T_pvar_handle handle) ;
int MPI_T_pvar_readreset(MPI_T_pvar_session session, MPI_T_pvar_handle handle, void *buf) ;
int MPI_T_category_get_num(int *num_cat) ;
int MPI_T_category_get_info(int cat_index, char *name, int *name_len, char *desc, int *desc_len,
                            int *num_cvars, int *num_pvars, int *num_categories) ;
int MPI_T_category_get_cvars(int cat_index, int len, int indices[]) ;
int MPI_T_category_get_pvars(int cat_index, int len, int indices[]) ;
int MPI_T_category_get_categories(int cat_index, int len, int indices[]) ;
int MPI_T_category_changed(int *stamp) ;
int MPI_T_cvar_get_index(const char *name, int *cvar_index) ;
int MPI_T_pvar_get_index(const char *name, int var_class, int *pvar_index) ;
int MPI_T_category_get_index(const char *name, int *cat_index) ;
 


 
 
int MPIX_Comm_failure_ack(MPI_Comm comm) ;
int MPIX_Comm_failure_get_acked(MPI_Comm comm, MPI_Group *failedgrp) ;
int MPIX_Comm_revoke(MPI_Comm comm) ;
int MPIX_Comm_shrink(MPI_Comm comm, MPI_Comm *newcomm) ;
int MPIX_Comm_agree(MPI_Comm comm, int *flag) ;


 


 
int PMPI_Send(const void *buf, int count, MPI_Datatype datatype, int dest, int tag,
              MPI_Comm comm)  ;
int PMPI_Recv(void *buf, int count, MPI_Datatype datatype, int source, int tag,
              MPI_Comm comm, MPI_Status *status)  ;
int PMPI_Get_count(const MPI_Status *status, MPI_Datatype datatype, int *count) ;
int PMPI_Bsend(const void *buf, int count, MPI_Datatype datatype, int dest, int tag,
               MPI_Comm comm)  ;
int PMPI_Ssend(const void *buf, int count, MPI_Datatype datatype, int dest, int tag,
               MPI_Comm comm)  ;
int PMPI_Rsend(const void *buf, int count, MPI_Datatype datatype, int dest, int tag,
               MPI_Comm comm)  ;
int PMPI_Buffer_attach(void *buffer, int size) ;
int PMPI_Buffer_detach(void *buffer_addr, int *size) ;
int PMPI_Isend(const void *buf, int count, MPI_Datatype datatype, int dest, int tag,
               MPI_Comm comm, MPI_Request *request)  ;
int PMPI_Ibsend(const void *buf, int count, MPI_Datatype datatype, int dest, int tag,
                MPI_Comm comm, MPI_Request *request)  ;
int PMPI_Issend(const void *buf, int count, MPI_Datatype datatype, int dest, int tag,
                MPI_Comm comm, MPI_Request *request)  ;
int PMPI_Irsend(const void *buf, int count, MPI_Datatype datatype, int dest, int tag,
                MPI_Comm comm, MPI_Request *request)  ;
int PMPI_Irecv(void *buf, int count, MPI_Datatype datatype, int source, int tag,
               MPI_Comm comm, MPI_Request *request)  ;
int PMPI_Wait(MPI_Request *request, MPI_Status *status) ;
int PMPI_Test(MPI_Request *request, int *flag, MPI_Status *status) ;
int PMPI_Request_free(MPI_Request *request) ;
int PMPI_Waitany(int count, MPI_Request array_of_requests[], int *indx, MPI_Status *status) ;
int PMPI_Testany(int count, MPI_Request array_of_requests[], int *indx, int *flag,
                 MPI_Status *status) ;
int PMPI_Waitall(int count, MPI_Request array_of_requests[], MPI_Status array_of_statuses[]) ;
int PMPI_Testall(int count, MPI_Request array_of_requests[], int *flag,
                 MPI_Status array_of_statuses[]) ;
int PMPI_Waitsome(int incount, MPI_Request array_of_requests[], int *outcount,
                  int array_of_indices[], MPI_Status array_of_statuses[]) ;
int PMPI_Testsome(int incount, MPI_Request array_of_requests[], int *outcount,
                  int array_of_indices[], MPI_Status array_of_statuses[]) ;
int PMPI_Iprobe(int source, int tag, MPI_Comm comm, int *flag, MPI_Status *status) ;
int PMPI_Probe(int source, int tag, MPI_Comm comm, MPI_Status *status) ;
int PMPI_Cancel(MPI_Request *request) ;
int PMPI_Test_cancelled(const MPI_Status *status, int *flag) ;
int PMPI_Send_init(const void *buf, int count, MPI_Datatype datatype, int dest, int tag,
                   MPI_Comm comm, MPI_Request *request)  ;
int PMPI_Bsend_init(const void *buf, int count, MPI_Datatype datatype, int dest, int tag,
                    MPI_Comm comm, MPI_Request *request)  ;
int PMPI_Ssend_init(const void *buf, int count, MPI_Datatype datatype, int dest, int tag,
                    MPI_Comm comm, MPI_Request *request)  ;
int PMPI_Rsend_init(const void *buf, int count, MPI_Datatype datatype, int dest, int tag,
                    MPI_Comm comm, MPI_Request *request)  ;
int PMPI_Recv_init(void *buf, int count, MPI_Datatype datatype, int source, int tag,
                   MPI_Comm comm, MPI_Request *request)  ;
int PMPI_Start(MPI_Request *request) ;
int PMPI_Startall(int count, MPI_Request array_of_requests[]) ;
int PMPI_Sendrecv(const void *sendbuf, int sendcount, MPI_Datatype sendtype, int dest,
                  int sendtag, void *recvbuf, int recvcount, MPI_Datatype recvtype,
                  int source, int recvtag, MPI_Comm comm, MPI_Status *status)
                  
                   ;
int PMPI_Sendrecv_replace(void *buf, int count, MPI_Datatype datatype, int dest,
                          int sendtag, int source, int recvtag, MPI_Comm comm,
                          MPI_Status *status)  ;
int PMPI_Type_contiguous(int count, MPI_Datatype oldtype, MPI_Datatype *newtype) ;
int PMPI_Type_vector(int count, int blocklength, int stride, MPI_Datatype oldtype,
                     MPI_Datatype *newtype) ;
int PMPI_Type_hvector(int count, int blocklength, MPI_Aint stride, MPI_Datatype oldtype,
                      MPI_Datatype *newtype) ;
int PMPI_Type_indexed(int count, const int *array_of_blocklengths,
                      const int *array_of_displacements, MPI_Datatype oldtype,
                      MPI_Datatype *newtype) ;
int PMPI_Type_hindexed(int count, int *array_of_blocklengths,
                       MPI_Aint *array_of_displacements, MPI_Datatype oldtype,
                       MPI_Datatype *newtype) ;
int PMPI_Type_struct(int count, int *array_of_blocklengths,
                     MPI_Aint *array_of_displacements,
                     MPI_Datatype *array_of_types, MPI_Datatype *newtype) ;
int PMPI_Address(void *location, MPI_Aint *address) ;
int PMPI_Type_extent(MPI_Datatype datatype, MPI_Aint *extent) ;
int PMPI_Type_size(MPI_Datatype datatype, int *size) ;
int PMPI_Type_lb(MPI_Datatype datatype, MPI_Aint *displacement) ;
int PMPI_Type_ub(MPI_Datatype datatype, MPI_Aint *displacement) ;
int PMPI_Type_commit(MPI_Datatype *datatype) ;
int PMPI_Type_free(MPI_Datatype *datatype) ;
int PMPI_Get_elements(const MPI_Status *status, MPI_Datatype datatype, int *count) ;
int PMPI_Pack(const void *inbuf, int incount, MPI_Datatype datatype, void *outbuf,
              int outsize, int *position, MPI_Comm comm)  ;
int PMPI_Unpack(const void *inbuf, int insize, int *position, void *outbuf, int outcount,
                MPI_Datatype datatype, MPI_Comm comm)  ;
int PMPI_Pack_size(int incount, MPI_Datatype datatype, MPI_Comm comm, int *size) ;
int PMPI_Barrier(MPI_Comm comm) ;
int PMPI_Bcast(void *buffer, int count, MPI_Datatype datatype, int root, MPI_Comm comm)
                ;
int PMPI_Gather(const void *sendbuf, int sendcount, MPI_Datatype sendtype, void *recvbuf,
                int recvcount, MPI_Datatype recvtype, int root, MPI_Comm comm)
                  ;
int PMPI_Gatherv(const void *sendbuf, int sendcount, MPI_Datatype sendtype, void *recvbuf,
                 const int *recvcounts, const int *displs, MPI_Datatype recvtype, int root,
                 MPI_Comm comm)
                   ;
int PMPI_Scatter(const void *sendbuf, int sendcount, MPI_Datatype sendtype, void *recvbuf,
                 int recvcount, MPI_Datatype recvtype, int root, MPI_Comm comm)
                   ;
int PMPI_Scatterv(const void *sendbuf, const int *sendcounts, const int *displs,
                  MPI_Datatype sendtype, void *recvbuf, int recvcount, MPI_Datatype recvtype,
                  int root, MPI_Comm comm)
                    ;
int PMPI_Allgather(const void *sendbuf, int sendcount, MPI_Datatype sendtype, void *recvbuf,
                   int recvcount, MPI_Datatype recvtype, MPI_Comm comm)
                     ;
int PMPI_Allgatherv(const void *sendbuf, int sendcount, MPI_Datatype sendtype, void *recvbuf,
                    const int *recvcounts, const int *displs, MPI_Datatype recvtype, MPI_Comm comm)
                      ;
int PMPI_Alltoall(const void *sendbuf, int sendcount, MPI_Datatype sendtype, void *recvbuf,
                  int recvcount, MPI_Datatype recvtype, MPI_Comm comm)
                    ;
int PMPI_Alltoallv(const void *sendbuf, const int *sendcounts, const int *sdispls,
                   MPI_Datatype sendtype, void *recvbuf, const int *recvcounts,
                   const int *rdispls, MPI_Datatype recvtype, MPI_Comm comm)
                     ;
int PMPI_Alltoallw(const void *sendbuf, const int sendcounts[], const int sdispls[],
                   const MPI_Datatype sendtypes[], void *recvbuf, const int recvcounts[],
                   const int rdispls[], const MPI_Datatype recvtypes[], MPI_Comm comm) ;
int PMPI_Exscan(const void *sendbuf, void *recvbuf, int count, MPI_Datatype datatype,
                MPI_Op op, MPI_Comm comm)
                  ;
int PMPI_Reduce(const void *sendbuf, void *recvbuf, int count, MPI_Datatype datatype,
                MPI_Op op, int root, MPI_Comm comm)
                  ;
int PMPI_Op_create(MPI_User_function *user_fn, int commute, MPI_Op *op) ;
int PMPI_Op_free(MPI_Op *op) ;
int PMPI_Allreduce(const void *sendbuf, void *recvbuf, int count, MPI_Datatype datatype,
                   MPI_Op op, MPI_Comm comm)
                     ;
int PMPI_Reduce_scatter(const void *sendbuf, void *recvbuf, const int recvcounts[],
                        MPI_Datatype datatype, MPI_Op op, MPI_Comm comm)
                          ;
int PMPI_Scan(const void *sendbuf, void *recvbuf, int count, MPI_Datatype datatype, MPI_Op op,
              MPI_Comm comm)
                ;
int PMPI_Group_size(MPI_Group group, int *size) ;
int PMPI_Group_rank(MPI_Group group, int *rank) ;
int PMPI_Group_translate_ranks(MPI_Group group1, int n, const int ranks1[], MPI_Group group2,
                               int ranks2[]) ;
int PMPI_Group_compare(MPI_Group group1, MPI_Group group2, int *result) ;
int PMPI_Comm_group(MPI_Comm comm, MPI_Group *group) ;
int PMPI_Group_union(MPI_Group group1, MPI_Group group2, MPI_Group *newgroup) ;
int PMPI_Group_intersection(MPI_Group group1, MPI_Group group2, MPI_Group *newgroup) ;
int PMPI_Group_difference(MPI_Group group1, MPI_Group group2, MPI_Group *newgroup) ;
int PMPI_Group_incl(MPI_Group group, int n, const int ranks[], MPI_Group *newgroup) ;
int PMPI_Group_excl(MPI_Group group, int n, const int ranks[], MPI_Group *newgroup) ;
int PMPI_Group_range_incl(MPI_Group group, int n, int ranges[][3], MPI_Group *newgroup) ;
int PMPI_Group_range_excl(MPI_Group group, int n, int ranges[][3], MPI_Group *newgroup) ;
int PMPI_Group_free(MPI_Group *group) ;
int PMPI_Comm_size(MPI_Comm comm, int *size) ;
int PMPI_Comm_rank(MPI_Comm comm, int *rank) ;
int PMPI_Comm_compare(MPI_Comm comm1, MPI_Comm comm2, int *result) ;
int PMPI_Comm_dup(MPI_Comm comm, MPI_Comm *newcomm) ;
int PMPI_Comm_dup_with_info(MPI_Comm comm, MPI_Info info, MPI_Comm *newcomm) ;
int PMPI_Comm_create(MPI_Comm comm, MPI_Group group, MPI_Comm *newcomm) ;
int PMPI_Comm_split(MPI_Comm comm, int color, int key, MPI_Comm *newcomm) ;
int PMPI_Comm_free(MPI_Comm *comm) ;
int PMPI_Comm_test_inter(MPI_Comm comm, int *flag) ;
int PMPI_Comm_remote_size(MPI_Comm comm, int *size) ;
int PMPI_Comm_remote_group(MPI_Comm comm, MPI_Group *group) ;
int PMPI_Intercomm_create(MPI_Comm local_comm, int local_leader, MPI_Comm peer_comm,
                          int remote_leader, int tag, MPI_Comm *newintercomm) ;
int PMPI_Intercomm_merge(MPI_Comm intercomm, int high, MPI_Comm *newintracomm) ;
int PMPI_Keyval_create(MPI_Copy_function *copy_fn, MPI_Delete_function *delete_fn,
                       int *keyval, void *extra_state) ;
int PMPI_Keyval_free(int *keyval) ;
int PMPI_Attr_put(MPI_Comm comm, int keyval, void *attribute_val) ;
int PMPI_Attr_get(MPI_Comm comm, int keyval, void *attribute_val, int *flag) ;
int PMPI_Attr_delete(MPI_Comm comm, int keyval) ;
int PMPI_Topo_test(MPI_Comm comm, int *status) ;
int PMPI_Cart_create(MPI_Comm comm_old, int ndims, const int dims[], const int periods[],
                     int reorder, MPI_Comm *comm_cart) ;
int PMPI_Dims_create(int nnodes, int ndims, int dims[]) ;
int PMPI_Graph_create(MPI_Comm comm_old, int nnodes, const int indx[], const int edges[],
                      int reorder, MPI_Comm *comm_graph) ;
int PMPI_Graphdims_get(MPI_Comm comm, int *nnodes, int *nedges) ;
int PMPI_Graph_get(MPI_Comm comm, int maxindex, int maxedges, int indx[], int edges[]) ;
int PMPI_Cartdim_get(MPI_Comm comm, int *ndims) ;
int PMPI_Cart_get(MPI_Comm comm, int maxdims, int dims[], int periods[], int coords[]) ;
int PMPI_Cart_rank(MPI_Comm comm, const int coords[], int *rank) ;
int PMPI_Cart_coords(MPI_Comm comm, int rank, int maxdims, int coords[]) ;
int PMPI_Graph_neighbors_count(MPI_Comm comm, int rank, int *nneighbors) ;
int PMPI_Graph_neighbors(MPI_Comm comm, int rank, int maxneighbors, int neighbors[]) ;
int PMPI_Cart_shift(MPI_Comm comm, int direction, int disp, int *rank_source, int *rank_dest) ;
int PMPI_Cart_sub(MPI_Comm comm, const int remain_dims[], MPI_Comm *newcomm) ;
int PMPI_Cart_map(MPI_Comm comm, int ndims, const int dims[], const int periods[], int *newrank) ;
int PMPI_Graph_map(MPI_Comm comm, int nnodes, const int indx[], const int edges[], int *newrank) ;
int PMPI_Get_processor_name(char *name, int *resultlen) ;
int PMPI_Get_version(int *version, int *subversion) ;
int PMPI_Get_library_version(char *version, int *resultlen) ;
int PMPI_Errhandler_create(MPI_Handler_function *function, MPI_Errhandler *errhandler) ;
int PMPI_Errhandler_set(MPI_Comm comm, MPI_Errhandler errhandler) ;
int PMPI_Errhandler_get(MPI_Comm comm, MPI_Errhandler *errhandler) ;
int PMPI_Errhandler_free(MPI_Errhandler *errhandler) ;
int PMPI_Error_string(int errorcode, char *string, int *resultlen) ;
int PMPI_Error_class(int errorcode, int *errorclass) ;
double PMPI_Wtime(void) ;
double PMPI_Wtick(void) ;
int PMPI_Init(int *argc, char ***argv) ;
int PMPI_Finalize(void) ;
int PMPI_Initialized(int *flag) ;
int PMPI_Abort(MPI_Comm comm, int errorcode) ;


 
int PMPI_Pcontrol(const int level, ...) ;

 
int PMPI_Close_port(const char *port_name) ;
int PMPI_Comm_accept(const char *port_name, MPI_Info info, int root, MPI_Comm comm,
                     MPI_Comm *newcomm) ;
int PMPI_Comm_connect(const char *port_name, MPI_Info info, int root, MPI_Comm comm,
                      MPI_Comm *newcomm) ;
int PMPI_Comm_disconnect(MPI_Comm *comm) ;
int PMPI_Comm_get_parent(MPI_Comm *parent) ;
int PMPI_Comm_join(int fd, MPI_Comm *intercomm) ;
int PMPI_Comm_spawn(const char *command, char *argv[], int maxprocs, MPI_Info info, int root,
                    MPI_Comm comm, MPI_Comm *intercomm, int array_of_errcodes[]) ;
int PMPI_Comm_spawn_multiple(int count, char *array_of_commands[], char **array_of_argv[],
                             const int array_of_maxprocs[], const MPI_Info array_of_info[],
                             int root, MPI_Comm comm, MPI_Comm *intercomm, int array_of_errcodes[]) ;
int PMPI_Lookup_name(const char *service_name, MPI_Info info, char *port_name) ;
int PMPI_Open_port(MPI_Info info, char *port_name) ;
int PMPI_Publish_name(const char *service_name, MPI_Info info, const char *port_name) ;
int PMPI_Unpublish_name(const char *service_name, MPI_Info info, const char *port_name) ;
int PMPI_Comm_set_info(MPI_Comm comm, MPI_Info info) ;
int PMPI_Comm_get_info(MPI_Comm comm, MPI_Info *info) ;

 
int PMPI_Accumulate(const void *origin_addr, int origin_count, MPI_Datatype origin_datatype,
                    int target_rank, MPI_Aint target_disp, int target_count,
                    MPI_Datatype target_datatype, MPI_Op op, MPI_Win win)
                     ;
int PMPI_Get(void *origin_addr, int origin_count, MPI_Datatype origin_datatype,
             int target_rank, MPI_Aint target_disp, int target_count,
             MPI_Datatype target_datatype, MPI_Win win)  ;
int PMPI_Put(const void *origin_addr, int origin_count, MPI_Datatype origin_datatype,
             int target_rank, MPI_Aint target_disp, int target_count,
             MPI_Datatype target_datatype, MPI_Win win)  ;
int PMPI_Win_complete(MPI_Win win) ;
int PMPI_Win_create(void *base, MPI_Aint size, int disp_unit, MPI_Info info, MPI_Comm comm,
                    MPI_Win *win) ;
int PMPI_Win_fence(int assert, MPI_Win win) ;
int PMPI_Win_free(MPI_Win *win) ;
int PMPI_Win_get_group(MPI_Win win, MPI_Group *group) ;
int PMPI_Win_lock(int lock_type, int rank, int assert, MPI_Win win) ;
int PMPI_Win_post(MPI_Group group, int assert, MPI_Win win) ;
int PMPI_Win_start(MPI_Group group, int assert, MPI_Win win) ;
int PMPI_Win_test(MPI_Win win, int *flag) ;
int PMPI_Win_unlock(int rank, MPI_Win win) ;
int PMPI_Win_wait(MPI_Win win) ;

 
int PMPI_Win_allocate(MPI_Aint size, int disp_unit, MPI_Info info, MPI_Comm comm, void *baseptr,
                      MPI_Win *win) ;
int PMPI_Win_allocate_shared(MPI_Aint size, int disp_unit, MPI_Info info, MPI_Comm comm,
                             void *baseptr, MPI_Win *win) ;
int PMPI_Win_shared_query(MPI_Win win, int rank, MPI_Aint *size, int *disp_unit, void *baseptr) ;
int PMPI_Win_create_dynamic(MPI_Info info, MPI_Comm comm, MPI_Win *win) ;
int PMPI_Win_attach(MPI_Win win, void *base, MPI_Aint size) ;
int PMPI_Win_detach(MPI_Win win, const void *base) ;
int PMPI_Win_get_info(MPI_Win win, MPI_Info *info_used) ;
int PMPI_Win_set_info(MPI_Win win, MPI_Info info) ;
int PMPI_Get_accumulate(const void *origin_addr, int origin_count,
                         MPI_Datatype origin_datatype, void *result_addr, int result_count,
                         MPI_Datatype result_datatype, int target_rank, MPI_Aint target_disp,
                         int target_count, MPI_Datatype target_datatype, MPI_Op op, MPI_Win win)
                         
                          ;
int PMPI_Fetch_and_op(const void *origin_addr, void *result_addr,
                       MPI_Datatype datatype, int target_rank, MPI_Aint target_disp,
                       MPI_Op op, MPI_Win win)
                        ;
int PMPI_Compare_and_swap(const void *origin_addr, const void *compare_addr,
                           void *result_addr, MPI_Datatype datatype, int target_rank,
                           MPI_Aint target_disp, MPI_Win win)
                           
                           
                            ;
int PMPI_Rput(const void *origin_addr, int origin_count,
               MPI_Datatype origin_datatype, int target_rank, MPI_Aint target_disp,
               int target_count, MPI_Datatype target_datatype, MPI_Win win,
               MPI_Request *request)
                ;
int PMPI_Rget(void *origin_addr, int origin_count,
               MPI_Datatype origin_datatype, int target_rank, MPI_Aint target_disp,
               int target_count, MPI_Datatype target_datatype, MPI_Win win,
               MPI_Request *request)
                ;
int PMPI_Raccumulate(const void *origin_addr, int origin_count,
                      MPI_Datatype origin_datatype, int target_rank, MPI_Aint target_disp,
                      int target_count, MPI_Datatype target_datatype, MPI_Op op, MPI_Win win,
                      MPI_Request *request)
                       ;
int PMPI_Rget_accumulate(const void *origin_addr, int origin_count,
                          MPI_Datatype origin_datatype, void *result_addr, int result_count,
                          MPI_Datatype result_datatype, int target_rank, MPI_Aint target_disp,
                          int target_count, MPI_Datatype target_datatype, MPI_Op op, MPI_Win win,
                          MPI_Request *request)
                          
                           ;
int PMPI_Win_lock_all(int assert, MPI_Win win) ;
int PMPI_Win_unlock_all(MPI_Win win) ;
int PMPI_Win_flush(int rank, MPI_Win win) ;
int PMPI_Win_flush_all(MPI_Win win) ;
int PMPI_Win_flush_local(int rank, MPI_Win win) ;
int PMPI_Win_flush_local_all(MPI_Win win) ;
int PMPI_Win_sync(MPI_Win win) ;
 
 
int PMPI_Add_error_class(int *errorclass) ;
int PMPI_Add_error_code(int errorclass, int *errorcode) ;
int PMPI_Add_error_string(int errorcode, const char *string) ;
int PMPI_Comm_call_errhandler(MPI_Comm comm, int errorcode) ;
int PMPI_Comm_create_keyval(MPI_Comm_copy_attr_function *comm_copy_attr_fn,
                            MPI_Comm_delete_attr_function *comm_delete_attr_fn, int *comm_keyval,
                            void *extra_state) ;
int PMPI_Comm_delete_attr(MPI_Comm comm, int comm_keyval) ;
int PMPI_Comm_free_keyval(int *comm_keyval) ;
int PMPI_Comm_get_attr(MPI_Comm comm, int comm_keyval, void *attribute_val, int *flag) ;
int PMPI_Comm_get_name(MPI_Comm comm, char *comm_name, int *resultlen) ;
int PMPI_Comm_set_attr(MPI_Comm comm, int comm_keyval, void *attribute_val) ;
int PMPI_Comm_set_name(MPI_Comm comm, const char *comm_name) ;
int PMPI_File_call_errhandler(MPI_File fh, int errorcode) ;
int PMPI_Grequest_complete(MPI_Request request) ;
int PMPI_Grequest_start(MPI_Grequest_query_function *query_fn, MPI_Grequest_free_function *free_fn,
                        MPI_Grequest_cancel_function *cancel_fn, void *extra_state,
                        MPI_Request *request) ;
int PMPI_Init_thread(int *argc, char ***argv, int required, int *provided) ;
int PMPI_Is_thread_main(int *flag) ;
int PMPI_Query_thread(int *provided) ;
int PMPI_Status_set_cancelled(MPI_Status *status, int flag) ;
int PMPI_Status_set_elements(MPI_Status *status, MPI_Datatype datatype, int count) ;
int PMPI_Type_create_keyval(MPI_Type_copy_attr_function *type_copy_attr_fn,
                            MPI_Type_delete_attr_function *type_delete_attr_fn,
                            int *type_keyval, void *extra_state) ;
int PMPI_Type_delete_attr(MPI_Datatype datatype, int type_keyval) ;
int PMPI_Type_dup(MPI_Datatype oldtype, MPI_Datatype *newtype) ;
int PMPI_Type_free_keyval(int *type_keyval) ;
int PMPI_Type_get_attr(MPI_Datatype datatype, int type_keyval, void *attribute_val, int *flag) ;
int PMPI_Type_get_contents(MPI_Datatype datatype, int max_integers, int max_addresses,
                           int max_datatypes, int array_of_integers[],
                           MPI_Aint array_of_addresses[], MPI_Datatype array_of_datatypes[]) ;
int PMPI_Type_get_envelope(MPI_Datatype datatype, int *num_integers, int *num_addresses,
                           int *num_datatypes, int *combiner) ;
int PMPI_Type_get_name(MPI_Datatype datatype, char *type_name, int *resultlen) ;
int PMPI_Type_set_attr(MPI_Datatype datatype, int type_keyval, void *attribute_val) ;
int PMPI_Type_set_name(MPI_Datatype datatype, const char *type_name) ;
int PMPI_Type_match_size(int typeclass, int size, MPI_Datatype *datatype) ;
int PMPI_Win_call_errhandler(MPI_Win win, int errorcode) ;
int PMPI_Win_create_keyval(MPI_Win_copy_attr_function *win_copy_attr_fn,
                           MPI_Win_delete_attr_function *win_delete_attr_fn, int *win_keyval,
                           void *extra_state) ;
int PMPI_Win_delete_attr(MPI_Win win, int win_keyval) ;
int PMPI_Win_free_keyval(int *win_keyval) ;
int PMPI_Win_get_attr(MPI_Win win, int win_keyval, void *attribute_val, int *flag) ;
int PMPI_Win_get_name(MPI_Win win, char *win_name, int *resultlen) ;
int PMPI_Win_set_attr(MPI_Win win, int win_keyval, void *attribute_val) ;
int PMPI_Win_set_name(MPI_Win win, const char *win_name) ;

int PMPI_Alloc_mem(MPI_Aint size, MPI_Info info, void *baseptr) ;
int PMPI_Comm_create_errhandler(MPI_Comm_errhandler_function *comm_errhandler_fn,
                                MPI_Errhandler *errhandler) ;
int PMPI_Comm_get_errhandler(MPI_Comm comm, MPI_Errhandler *errhandler) ;
int PMPI_Comm_set_errhandler(MPI_Comm comm, MPI_Errhandler errhandler) ;
int PMPI_File_create_errhandler(MPI_File_errhandler_function *file_errhandler_fn,
                                MPI_Errhandler *errhandler) ;
int PMPI_File_get_errhandler(MPI_File file, MPI_Errhandler *errhandler) ;
int PMPI_File_set_errhandler(MPI_File file, MPI_Errhandler errhandler) ;
int PMPI_Finalized(int *flag) ;
int PMPI_Free_mem(void *base) ;
int PMPI_Get_address(const void *location, MPI_Aint *address) ;
int PMPI_Info_create(MPI_Info *info) ;
int PMPI_Info_delete(MPI_Info info, const char *key) ;
int PMPI_Info_dup(MPI_Info info, MPI_Info *newinfo) ;
int PMPI_Info_free(MPI_Info *info) ;
int PMPI_Info_get(MPI_Info info, const char *key, int valuelen, char *value, int *flag) ;
int PMPI_Info_get_nkeys(MPI_Info info, int *nkeys) ;
int PMPI_Info_get_nthkey(MPI_Info info, int n, char *key) ;
int PMPI_Info_get_valuelen(MPI_Info info, const char *key, int *valuelen, int *flag) ;
int PMPI_Info_set(MPI_Info info, const char *key, const char *value) ;
int PMPI_Pack_external(const char datarep[], const void *inbuf, int incount,
                       MPI_Datatype datatype, void *outbuf, MPI_Aint outsize, MPI_Aint *position)
                        ;
int PMPI_Pack_external_size(const char datarep[], int incount, MPI_Datatype datatype,
                            MPI_Aint *size) ;
int PMPI_Request_get_status(MPI_Request request, int *flag, MPI_Status *status) ;
int PMPI_Status_c2f(const MPI_Status *c_status, MPI_Fint *f_status) ;
int PMPI_Status_f2c(const MPI_Fint *f_status, MPI_Status *c_status) ;
int PMPI_Type_create_darray(int size, int rank, int ndims, const int array_of_gsizes[],
                            const int array_of_distribs[], const int array_of_dargs[],
                            const int array_of_psizes[], int order, MPI_Datatype oldtype,
                            MPI_Datatype *newtype) ;
int PMPI_Type_create_hindexed(int count, const int array_of_blocklengths[],
                              const MPI_Aint array_of_displacements[], MPI_Datatype oldtype,
                              MPI_Datatype *newtype) ;
int PMPI_Type_create_hvector(int count, int blocklength, MPI_Aint stride, MPI_Datatype oldtype,
                             MPI_Datatype *newtype) ;
int PMPI_Type_create_indexed_block(int count, int blocklength, const int array_of_displacements[],
                                   MPI_Datatype oldtype, MPI_Datatype *newtype) ;
int PMPI_Type_create_hindexed_block(int count, int blocklength,
                                    const MPI_Aint array_of_displacements[],
                                    MPI_Datatype oldtype, MPI_Datatype *newtype) ;
int PMPI_Type_create_resized(MPI_Datatype oldtype, MPI_Aint lb, MPI_Aint extent,
                             MPI_Datatype *newtype) ;
int PMPI_Type_create_struct(int count, const int array_of_blocklengths[],
                            const MPI_Aint array_of_displacements[],
                            const MPI_Datatype array_of_types[], MPI_Datatype *newtype) ;
int PMPI_Type_create_subarray(int ndims, const int array_of_sizes[],
                              const int array_of_subsizes[], const int array_of_starts[],
                              int order, MPI_Datatype oldtype, MPI_Datatype *newtype) ;
int PMPI_Type_get_extent(MPI_Datatype datatype, MPI_Aint *lb, MPI_Aint *extent) ;
int PMPI_Type_get_true_extent(MPI_Datatype datatype, MPI_Aint *true_lb, MPI_Aint *true_extent) ;
int PMPI_Unpack_external(const char datarep[], const void *inbuf, MPI_Aint insize,
                         MPI_Aint *position, void *outbuf, int outcount, MPI_Datatype datatype)
                          ;
int PMPI_Win_create_errhandler(MPI_Win_errhandler_function *win_errhandler_fn,
                               MPI_Errhandler *errhandler) ;
int PMPI_Win_get_errhandler(MPI_Win win, MPI_Errhandler *errhandler) ;
int PMPI_Win_set_errhandler(MPI_Win win, MPI_Errhandler errhandler) ;



 
int PMPI_Type_create_f90_integer(int r, MPI_Datatype *newtype) ;
int PMPI_Type_create_f90_real(int p, int r, MPI_Datatype *newtype) ;
int PMPI_Type_create_f90_complex(int p, int r, MPI_Datatype *newtype) ;

int PMPI_Reduce_local(const void *inbuf, void *inoutbuf, int count, MPI_Datatype datatype,
                      MPI_Op op)
                        ;
int PMPI_Op_commutative(MPI_Op op, int *commute) ;
int PMPI_Reduce_scatter_block(const void *sendbuf, void *recvbuf, int recvcount,
                              MPI_Datatype datatype, MPI_Op op, MPI_Comm comm)
                              
                               ;
int PMPI_Dist_graph_create_adjacent(MPI_Comm comm_old, int indegree, const int sources[],
                                    const int sourceweights[], int outdegree,
                                    const int destinations[], const int destweights[],
                                    MPI_Info info, int reorder, MPI_Comm *comm_dist_graph) ;
int PMPI_Dist_graph_create(MPI_Comm comm_old, int n, const int sources[], const int degrees[],
                           const int destinations[], const int weights[], MPI_Info info,
                           int reorder, MPI_Comm *comm_dist_graph) ;
int PMPI_Dist_graph_neighbors_count(MPI_Comm comm, int *indegree, int *outdegree, int *weighted) ;
int PMPI_Dist_graph_neighbors(MPI_Comm comm, int maxindegree, int sources[], int sourceweights[],
                              int maxoutdegree, int destinations[], int destweights[]) ;

 
int PMPI_Improbe(int source, int tag, MPI_Comm comm, int *flag, MPI_Message *message,
                 MPI_Status *status) ;
int PMPI_Imrecv(void *buf, int count, MPI_Datatype datatype, MPI_Message *message,
                MPI_Request *request)  ;
int PMPI_Mprobe(int source, int tag, MPI_Comm comm, MPI_Message *message, MPI_Status *status) ;
int PMPI_Mrecv(void *buf, int count, MPI_Datatype datatype, MPI_Message *message,
               MPI_Status *status)  ;

 
int PMPI_Comm_idup(MPI_Comm comm, MPI_Comm *newcomm, MPI_Request *request) ;
int PMPI_Ibarrier(MPI_Comm comm, MPI_Request *request) ;
int PMPI_Ibcast(void *buffer, int count, MPI_Datatype datatype, int root, MPI_Comm comm,
                MPI_Request *request)  ;
int PMPI_Igather(const void *sendbuf, int sendcount, MPI_Datatype sendtype, void *recvbuf,
                 int recvcount, MPI_Datatype recvtype, int root, MPI_Comm comm,
                 MPI_Request *request)
                   ;
int PMPI_Igatherv(const void *sendbuf, int sendcount, MPI_Datatype sendtype, void *recvbuf,
                  const int recvcounts[], const int displs[], MPI_Datatype recvtype, int root,
                  MPI_Comm comm, MPI_Request *request)
                    ;
int PMPI_Iscatter(const void *sendbuf, int sendcount, MPI_Datatype sendtype, void *recvbuf,
                  int recvcount, MPI_Datatype recvtype, int root, MPI_Comm comm,
                  MPI_Request *request)
                    ;
int PMPI_Iscatterv(const void *sendbuf, const int sendcounts[], const int displs[],
                   MPI_Datatype sendtype, void *recvbuf, int recvcount, MPI_Datatype recvtype,
                   int root, MPI_Comm comm, MPI_Request *request)
                     ;
int PMPI_Iallgather(const void *sendbuf, int sendcount, MPI_Datatype sendtype, void *recvbuf,
                    int recvcount, MPI_Datatype recvtype, MPI_Comm comm, MPI_Request *request)
                      ;
int PMPI_Iallgatherv(const void *sendbuf, int sendcount, MPI_Datatype sendtype, void *recvbuf,
                     const int recvcounts[], const int displs[], MPI_Datatype recvtype,
                     MPI_Comm comm, MPI_Request *request)
                       ;
int PMPI_Ialltoall(const void *sendbuf, int sendcount, MPI_Datatype sendtype, void *recvbuf,
                   int recvcount, MPI_Datatype recvtype, MPI_Comm comm, MPI_Request *request)
                     ;
int PMPI_Ialltoallv(const void *sendbuf, const int sendcounts[], const int sdispls[],
                    MPI_Datatype sendtype, void *recvbuf, const int recvcounts[],
                    const int rdispls[], MPI_Datatype recvtype, MPI_Comm comm,
                    MPI_Request *request)
                      ;
int PMPI_Ialltoallw(const void *sendbuf, const int sendcounts[], const int sdispls[],
                    const MPI_Datatype sendtypes[], void *recvbuf, const int recvcounts[],
                    const int rdispls[], const MPI_Datatype recvtypes[], MPI_Comm comm,
                    MPI_Request *request) ;
int PMPI_Ireduce(const void *sendbuf, void *recvbuf, int count, MPI_Datatype datatype,
                 MPI_Op op, int root, MPI_Comm comm, MPI_Request *request)
                   ;
int PMPI_Iallreduce(const void *sendbuf, void *recvbuf, int count, MPI_Datatype datatype,
                    MPI_Op op, MPI_Comm comm, MPI_Request *request)
                      ;
int PMPI_Ireduce_scatter(const void *sendbuf, void *recvbuf, const int recvcounts[],
                         MPI_Datatype datatype, MPI_Op op, MPI_Comm comm, MPI_Request *request)
                         
                          ;
int PMPI_Ireduce_scatter_block(const void *sendbuf, void *recvbuf, int recvcount,
                               MPI_Datatype datatype, MPI_Op op, MPI_Comm comm,
                               MPI_Request *request)
                               
                                ;
int PMPI_Iscan(const void *sendbuf, void *recvbuf, int count, MPI_Datatype datatype, MPI_Op op,
               MPI_Comm comm, MPI_Request *request)
                 ;
int PMPI_Iexscan(const void *sendbuf, void *recvbuf, int count, MPI_Datatype datatype,
                 MPI_Op op, MPI_Comm comm, MPI_Request *request)
                   ;

 
int PMPI_Ineighbor_allgather(const void *sendbuf, int sendcount, MPI_Datatype sendtype,
                             void *recvbuf, int recvcount, MPI_Datatype recvtype,
                             MPI_Comm comm, MPI_Request *request)
                             
                              ;
int PMPI_Ineighbor_allgatherv(const void *sendbuf, int sendcount, MPI_Datatype sendtype,
                              void *recvbuf, const int recvcounts[], const int displs[],
                              MPI_Datatype recvtype, MPI_Comm comm, MPI_Request *request)
                              
                               ;
int PMPI_Ineighbor_alltoall(const void *sendbuf, int sendcount, MPI_Datatype sendtype,
                            void *recvbuf, int recvcount, MPI_Datatype recvtype, MPI_Comm comm,
                            MPI_Request *request)
                            
                             ;
int PMPI_Ineighbor_alltoallv(const void *sendbuf, const int sendcounts[], const int sdispls[],
                             MPI_Datatype sendtype, void *recvbuf, const int recvcounts[],
                             const int rdispls[], MPI_Datatype recvtype, MPI_Comm comm,
                             MPI_Request *request)
                             
                              ;
int PMPI_Ineighbor_alltoallw(const void *sendbuf, const int sendcounts[],
                             const MPI_Aint sdispls[], const MPI_Datatype sendtypes[],
                             void *recvbuf, const int recvcounts[], const MPI_Aint rdispls[],
                             const MPI_Datatype recvtypes[], MPI_Comm comm, MPI_Request *request) ;
int PMPI_Neighbor_allgather(const void *sendbuf, int sendcount, MPI_Datatype sendtype,
                            void *recvbuf, int recvcount, MPI_Datatype recvtype, MPI_Comm comm)
                            
                             ;
int PMPI_Neighbor_allgatherv(const void *sendbuf, int sendcount, MPI_Datatype sendtype,
                             void *recvbuf, const int recvcounts[], const int displs[],
                             MPI_Datatype recvtype, MPI_Comm comm)
                             
                              ;
int PMPI_Neighbor_alltoall(const void *sendbuf, int sendcount, MPI_Datatype sendtype,
                           void *recvbuf, int recvcount, MPI_Datatype recvtype, MPI_Comm comm)
                           
                            ;
int PMPI_Neighbor_alltoallv(const void *sendbuf, const int sendcounts[], const int sdispls[],
                            MPI_Datatype sendtype, void *recvbuf, const int recvcounts[],
                            const int rdispls[], MPI_Datatype recvtype, MPI_Comm comm)
                            
                             ;
int PMPI_Neighbor_alltoallw(const void *sendbuf, const int sendcounts[], const MPI_Aint sdispls[],
                            const MPI_Datatype sendtypes[], void *recvbuf, const int recvcounts[],
                            const MPI_Aint rdispls[], const MPI_Datatype recvtypes[],
                            MPI_Comm comm) ;

 
int PMPI_Comm_split_type(MPI_Comm comm, int split_type, int key, MPI_Info info, MPI_Comm *newcomm) ;

 
int PMPI_Comm_create_group(MPI_Comm comm, MPI_Group group, int tag, MPI_Comm *newcomm) ;

 
int PMPI_Get_elements_x(const MPI_Status *status, MPI_Datatype datatype, MPI_Count *count) ;
int PMPI_Status_set_elements_x(MPI_Status *status, MPI_Datatype datatype, MPI_Count count) ;
int PMPI_Type_get_extent_x(MPI_Datatype datatype, MPI_Count *lb, MPI_Count *extent) ;
int PMPI_Type_get_true_extent_x(MPI_Datatype datatype, MPI_Count *lb, MPI_Count *extent) ;
int PMPI_Type_size_x(MPI_Datatype datatype, MPI_Count *size) ;

 
MPI_Aint PMPI_Aint_add(MPI_Aint base, MPI_Aint disp) ;
MPI_Aint PMPI_Aint_diff(MPI_Aint addr1, MPI_Aint addr2) ;

 

 
 
int PMPI_T_init_thread(int required, int *provided) ;
int PMPI_T_finalize(void) ;
int PMPI_T_enum_get_info(MPI_T_enum enumtype, int *num, char *name, int *name_len) ;
int PMPI_T_enum_get_item(MPI_T_enum enumtype, int indx, int *value, char *name, int *name_len) ;
int PMPI_T_cvar_get_num(int *num_cvar) ;
int PMPI_T_cvar_get_info(int cvar_index, char *name, int *name_len, int *verbosity,
                         MPI_Datatype *datatype, MPI_T_enum *enumtype, char *desc, int *desc_len,
                         int *binding, int *scope) ;
int PMPI_T_cvar_handle_alloc(int cvar_index, void *obj_handle, MPI_T_cvar_handle *handle,
                             int *count) ;
int PMPI_T_cvar_handle_free(MPI_T_cvar_handle *handle) ;
int PMPI_T_cvar_read(MPI_T_cvar_handle handle, void *buf) ;
int PMPI_T_cvar_write(MPI_T_cvar_handle handle, const void *buf) ;
int PMPI_T_pvar_get_num(int *num_pvar) ;
int PMPI_T_pvar_get_info(int pvar_index, char *name, int *name_len, int *verbosity, int *var_class,
                         MPI_Datatype *datatype, MPI_T_enum *enumtype, char *desc, int *desc_len,
                         int *binding, int *readonly, int *continuous, int *atomic) ;
int PMPI_T_pvar_session_create(MPI_T_pvar_session *session) ;
int PMPI_T_pvar_session_free(MPI_T_pvar_session *session) ;
int PMPI_T_pvar_handle_alloc(MPI_T_pvar_session session, int pvar_index, void *obj_handle,
                             MPI_T_pvar_handle *handle, int *count) ;
int PMPI_T_pvar_handle_free(MPI_T_pvar_session session, MPI_T_pvar_handle *handle) ;
int PMPI_T_pvar_start(MPI_T_pvar_session session, MPI_T_pvar_handle handle) ;
int PMPI_T_pvar_stop(MPI_T_pvar_session session, MPI_T_pvar_handle handle) ;
int PMPI_T_pvar_read(MPI_T_pvar_session session, MPI_T_pvar_handle handle, void *buf) ;
int PMPI_T_pvar_write(MPI_T_pvar_session session, MPI_T_pvar_handle handle, const void *buf) ;
int PMPI_T_pvar_reset(MPI_T_pvar_session session, MPI_T_pvar_handle handle) ;
int PMPI_T_pvar_readreset(MPI_T_pvar_session session, MPI_T_pvar_handle handle, void *buf) ;
int PMPI_T_category_get_num(int *num_cat) ;
int PMPI_T_category_get_info(int cat_index, char *name, int *name_len, char *desc, int *desc_len,
                             int *num_cvars, int *num_pvars, int *num_categories) ;
int PMPI_T_category_get_cvars(int cat_index, int len, int indices[]) ;
int PMPI_T_category_get_pvars(int cat_index, int len, int indices[]) ;
int PMPI_T_category_get_categories(int cat_index, int len, int indices[]) ;
int PMPI_T_category_changed(int *stamp) ;
int PMPI_T_cvar_get_index(const char *name, int *cvar_index) ;
int PMPI_T_pvar_get_index(const char *name, int var_class, int *pvar_index) ;
int PMPI_T_category_get_index(const char *name, int *cat_index) ;
 


 
 
int PMPIX_Comm_failure_ack(MPI_Comm comm) ;
int PMPIX_Comm_failure_get_acked(MPI_Comm comm, MPI_Group *failedgrp) ;
int PMPIX_Comm_revoke(MPI_Comm comm) ;
int PMPIX_Comm_shrink(MPI_Comm comm, MPI_Comm *newcomm) ;
int PMPIX_Comm_agree(MPI_Comm comm, int *flag) ;

 

 













 
 











































 

 














 
 










































 
 




 


 
 
 


 

 

 
 






 

 
 


 


 


 
 
int MPI_File_open(MPI_Comm comm, const char *filename, int amode, MPI_Info info, MPI_File *fh) ;
int MPI_File_close(MPI_File *fh) ;
int MPI_File_delete(const char *filename, MPI_Info info) ;
int MPI_File_set_size(MPI_File fh, MPI_Offset size) ;
int MPI_File_preallocate(MPI_File fh, MPI_Offset size) ;
int MPI_File_get_size(MPI_File fh, MPI_Offset *size) ;
int MPI_File_get_group(MPI_File fh, MPI_Group *group) ;
int MPI_File_get_amode(MPI_File fh, int *amode) ;
int MPI_File_set_info(MPI_File fh, MPI_Info info) ;
int MPI_File_get_info(MPI_File fh, MPI_Info *info_used) ;

 
int MPI_File_set_view(MPI_File fh, MPI_Offset disp, MPI_Datatype etype, MPI_Datatype filetype,
                      const char *datarep, MPI_Info info) ;
int MPI_File_get_view(MPI_File fh, MPI_Offset *disp, MPI_Datatype *etype, MPI_Datatype *filetype,
                      char *datarep) ;

 
int MPI_File_read_at(MPI_File fh, MPI_Offset offset, void *buf, int count, MPI_Datatype datatype,
                     MPI_Status *status)  ;
int MPI_File_read_at_all(MPI_File fh, MPI_Offset offset, void * buf, int count,
                         MPI_Datatype datatype, MPI_Status *status)
     ;
int MPI_File_write_at(MPI_File fh, MPI_Offset offset, const void * buf, int count,
                      MPI_Datatype datatype, MPI_Status *status)
     ;
int MPI_File_write_at_all(MPI_File fh, MPI_Offset offset, const void *buf, int count,
                          MPI_Datatype datatype, MPI_Status *status)
     ;



  
int MPI_File_iread_at(MPI_File fh, MPI_Offset offset, void *buf, int count, MPI_Datatype datatype,
                      MPI_Request *request)  ;
int MPI_File_iwrite_at(MPI_File fh, MPI_Offset offset, const void *buf, int count,
                       MPI_Datatype datatype, MPI_Request *request)
     ;

 
int MPI_File_read(MPI_File fh, void *buf, int count, MPI_Datatype datatype, MPI_Status *status)
     ;
int MPI_File_read_all(MPI_File fh, void *buf, int count, MPI_Datatype datatype, MPI_Status *status)
     ;
int MPI_File_write(MPI_File fh, const void *buf, int count, MPI_Datatype datatype,
                   MPI_Status *status)  ;
int MPI_File_write_all(MPI_File fh, const void *buf, int count, MPI_Datatype datatype,
                       MPI_Status *status)  ;



  

int MPI_File_iread(MPI_File fh, void *buf, int count, MPI_Datatype datatype, MPI_Request *request)
     ;
int MPI_File_iwrite(MPI_File fh, const void *buf, int count, MPI_Datatype datatype,
                    MPI_Request *request)  ;

int MPI_File_seek(MPI_File fh, MPI_Offset offset, int whence) ;
int MPI_File_get_position(MPI_File fh, MPI_Offset *offset) ;
int MPI_File_get_byte_offset(MPI_File fh, MPI_Offset offset, MPI_Offset *disp) ;

 
int MPI_File_read_shared(MPI_File fh, void *buf, int count, MPI_Datatype datatype,
                         MPI_Status *status)  ;
int MPI_File_write_shared(MPI_File fh, const void *buf, int count, MPI_Datatype datatype,
                          MPI_Status *status)  ;
int MPI_File_iread_shared(MPI_File fh, void *buf, int count, MPI_Datatype datatype,
                          MPI_Request *request)  ;
int MPI_File_iwrite_shared(MPI_File fh, const void *buf, int count, MPI_Datatype datatype,
                           MPI_Request *request)  ;
int MPI_File_read_ordered(MPI_File fh, void *buf, int count, MPI_Datatype datatype,
                          MPI_Status *status)  ;
int MPI_File_write_ordered(MPI_File fh, const void *buf, int count, MPI_Datatype datatype,
                           MPI_Status *status)  ;
int MPI_File_seek_shared(MPI_File fh, MPI_Offset offset, int whence) ;
int MPI_File_get_position_shared(MPI_File fh, MPI_Offset *offset) ;

 
int MPI_File_read_at_all_begin(MPI_File fh, MPI_Offset offset, void *buf, int count,
                               MPI_Datatype datatype)  ;
int MPI_File_read_at_all_end(MPI_File fh, void *buf, MPI_Status *status) ;
int MPI_File_write_at_all_begin(MPI_File fh, MPI_Offset offset, const void *buf, int count,
                                MPI_Datatype datatype)  ;
int MPI_File_write_at_all_end(MPI_File fh, const void *buf, MPI_Status *status) ;
int MPI_File_read_all_begin(MPI_File fh, void *buf, int count, MPI_Datatype datatype)
     ;
int MPI_File_read_all_end(MPI_File fh, void *buf, MPI_Status *status) ;
int MPI_File_write_all_begin(MPI_File fh, const void *buf, int count, MPI_Datatype datatype)
     ;
int MPI_File_write_all_end(MPI_File fh, const void *buf, MPI_Status *status) ;
int MPI_File_read_ordered_begin(MPI_File fh, void *buf, int count, MPI_Datatype datatype)
     ;
int MPI_File_read_ordered_end(MPI_File fh, void *buf, MPI_Status *status) ;
int MPI_File_write_ordered_begin(MPI_File fh, const void *buf, int count, MPI_Datatype datatype)
     ;
int MPI_File_write_ordered_end(MPI_File fh, const void *buf, MPI_Status *status) ;

 
int MPI_File_get_type_extent(MPI_File fh, MPI_Datatype datatype, MPI_Aint *extent) ;

 
int MPI_Register_datarep(const char *datarep, MPI_Datarep_conversion_function *read_conversion_fn,
			 MPI_Datarep_conversion_function *write_conversion_fn,
			 MPI_Datarep_extent_function *dtype_file_extent_fn, void *extra_state) ;

 
int MPI_File_set_atomicity(MPI_File fh, int flag) ;
int MPI_File_get_atomicity(MPI_File fh, int *flag) ;
int MPI_File_sync(MPI_File fh) ;

 

 
int MPI_File_iread_at_all(MPI_File fh, MPI_Offset offset, void *buf, int count,
                           MPI_Datatype datatype, MPI_Request *request)
     ;
int MPI_File_iwrite_at_all(MPI_File fh, MPI_Offset offset, const void *buf, int count,
                            MPI_Datatype datatype, MPI_Request *request)
     ;
int MPI_File_iread_all(MPI_File fh, void *buf, int count, MPI_Datatype datatype,
                        MPI_Request *request)
     ;
int MPI_File_iwrite_all(MPI_File fh, const void *buf, int count, MPI_Datatype datatype,
                         MPI_Request *request)
     ;
 




 
 
 
MPI_File MPI_File_f2c(MPI_Fint file) ;
MPI_Fint MPI_File_c2f(MPI_File file) ;



 



 


 
int PMPI_File_open(MPI_Comm, const char *, int, MPI_Info, MPI_File *) ;
int PMPI_File_close(MPI_File *) ;
int PMPI_File_delete(const char *, MPI_Info) ;
int PMPI_File_set_size(MPI_File, MPI_Offset) ;
int PMPI_File_preallocate(MPI_File, MPI_Offset) ;
int PMPI_File_get_size(MPI_File, MPI_Offset *) ;
int PMPI_File_get_group(MPI_File, MPI_Group *) ;
int PMPI_File_get_amode(MPI_File, int *) ;
int PMPI_File_set_info(MPI_File, MPI_Info) ;
int PMPI_File_get_info(MPI_File, MPI_Info *) ;

 
int PMPI_File_set_view(MPI_File, MPI_Offset, 
    MPI_Datatype, MPI_Datatype, const char *, MPI_Info) ;
int PMPI_File_get_view(MPI_File, MPI_Offset *, 
      MPI_Datatype *, MPI_Datatype *, char *) ;

 
int PMPI_File_read_at(MPI_File, MPI_Offset, void *,
	      int, MPI_Datatype, MPI_Status *)
               ;
int PMPI_File_read_at_all(MPI_File, MPI_Offset, void *,
	      int, MPI_Datatype, MPI_Status *)
               ;
int PMPI_File_write_at(MPI_File, MPI_Offset, const void *,
	      int, MPI_Datatype, MPI_Status *)
               ;
int PMPI_File_write_at_all(MPI_File, MPI_Offset, const void *,
	      int, MPI_Datatype, MPI_Status *)
               ;



  

int PMPI_File_iread_at(MPI_File, MPI_Offset, void *,
	      int, MPI_Datatype, MPI_Request *)
               ;
int PMPI_File_iwrite_at(MPI_File, MPI_Offset, const void *,
	      int, MPI_Datatype, MPI_Request *)
               ;

 
int PMPI_File_read(MPI_File, void *, int, MPI_Datatype, MPI_Status *)
                    ;
int PMPI_File_read_all(MPI_File, void *, int, MPI_Datatype, MPI_Status *)
                        ;
int PMPI_File_write(MPI_File, const void *, int, MPI_Datatype, MPI_Status *)
                     ;
int PMPI_File_write_all(MPI_File, const void *, int, MPI_Datatype, MPI_Status *)
                         ;



  

int PMPI_File_iread(MPI_File, void *, int, MPI_Datatype, MPI_Request *)
                     ;
int PMPI_File_iwrite(MPI_File, const void *, int, MPI_Datatype, MPI_Request *)
                      ;

int PMPI_File_seek(MPI_File, MPI_Offset, int) ;
int PMPI_File_get_position(MPI_File, MPI_Offset *) ;
int PMPI_File_get_byte_offset(MPI_File, MPI_Offset, MPI_Offset *) ;

 
int PMPI_File_read_shared(MPI_File, void *, int, MPI_Datatype, MPI_Status *)
                           ;
int PMPI_File_write_shared(MPI_File, const void *, int, MPI_Datatype, MPI_Status *)
                            ;
int PMPI_File_iread_shared(MPI_File, void *, int, 
			   MPI_Datatype, MPI_Request *)
                            ;
int PMPI_File_iwrite_shared(MPI_File, const void *, int,
			    MPI_Datatype, MPI_Request *)
                             ;
int PMPI_File_read_ordered(MPI_File, void *, int, MPI_Datatype, MPI_Status *)
                            ;
int PMPI_File_write_ordered(MPI_File, const void *, int, MPI_Datatype, MPI_Status *)
                             ;
int PMPI_File_seek_shared(MPI_File, MPI_Offset, int) ;
int PMPI_File_get_position_shared(MPI_File, MPI_Offset *) ;

 
int PMPI_File_read_at_all_begin(MPI_File, MPI_Offset, void *,
                               int, MPI_Datatype)
                                ;
int PMPI_File_read_at_all_end(MPI_File, void *, MPI_Status *) ;
int PMPI_File_write_at_all_begin(MPI_File, MPI_Offset, const void *,
                                 int, MPI_Datatype)
                                  ;
int PMPI_File_write_at_all_end(MPI_File, const void *, MPI_Status *) ;
int PMPI_File_read_all_begin(MPI_File, void *, int, MPI_Datatype)
                              ;
int PMPI_File_read_all_end(MPI_File, void *, MPI_Status *) ;
int PMPI_File_write_all_begin(MPI_File, const void *, int, MPI_Datatype)
                               ;
int PMPI_File_write_all_end(MPI_File, const void *, MPI_Status *) ;
int PMPI_File_read_ordered_begin(MPI_File, void *, int, MPI_Datatype)
                                  ;
int PMPI_File_read_ordered_end(MPI_File, void *, MPI_Status *) ;
int PMPI_File_write_ordered_begin(MPI_File, const void *, int, MPI_Datatype)
                                   ;
int PMPI_File_write_ordered_end(MPI_File, const void *, MPI_Status *) ;

 
int PMPI_File_get_type_extent(MPI_File, MPI_Datatype, MPI_Aint *) ;

 
int PMPI_Register_datarep(const char *,
			 MPI_Datarep_conversion_function *,
			 MPI_Datarep_conversion_function *,
			 MPI_Datarep_extent_function *,
			 void *) ;

 
int PMPI_File_set_atomicity(MPI_File, int) ;
int PMPI_File_get_atomicity(MPI_File, int *) ;
int PMPI_File_sync(MPI_File) ;

 

 
int PMPI_File_iread_at_all(MPI_File fh, MPI_Offset offset, void *buf, int count,
                            MPI_Datatype datatype, MPI_Request *request)
     ;
int PMPI_File_iwrite_at_all(MPI_File fh, MPI_Offset offset, const void *buf, int count,
                             MPI_Datatype datatype, MPI_Request *request)
     ;
int PMPI_File_iread_all(MPI_File fh, void *buf, int count, MPI_Datatype datatype,
                         MPI_Request *request)
     ;
int PMPI_File_iwrite_all(MPI_File fh, const void *buf, int count, MPI_Datatype datatype,
                          MPI_Request *request)
     ;


 
MPI_File PMPI_File_f2c(MPI_Fint) ;
MPI_Fint PMPI_File_c2f(MPI_File) ;


 





 
typedef int MPIX_Grequest_class;
int MPIX_Grequest_class_create(MPI_Grequest_query_function *query_fn,
                               MPI_Grequest_free_function *free_fn,
                               MPI_Grequest_cancel_function *cancel_fn,
                               MPIX_Grequest_poll_function *poll_fn,
                               MPIX_Grequest_wait_function *wait_fn,
                               MPIX_Grequest_class *greq_class) ;
int MPIX_Grequest_class_allocate(MPIX_Grequest_class greq_class, void *extra_state,
                                 MPI_Request *request) ;
int MPIX_Grequest_start(MPI_Grequest_query_function *query_fn,
                        MPI_Grequest_free_function *free_fn,
                        MPI_Grequest_cancel_function *cancel_fn,
                        MPIX_Grequest_poll_function *poll_fn,
                        MPIX_Grequest_wait_function *wait_fn, void *extra_state,
                        MPI_Request *request) ;

 
int PMPIX_Grequest_class_create(MPI_Grequest_query_function *query_fn,
                                MPI_Grequest_free_function *free_fn,
                                MPI_Grequest_cancel_function *cancel_fn,
                                MPIX_Grequest_poll_function *poll_fn,
                                MPIX_Grequest_wait_function *wait_fn,
                                MPIX_Grequest_class *greq_class) ;
int PMPIX_Grequest_class_allocate(MPIX_Grequest_class greq_class, void *extra_state,
                                  MPI_Request *request) ;
int PMPIX_Grequest_start(MPI_Grequest_query_function *query_fn,
                         MPI_Grequest_free_function *free_fn,
                         MPI_Grequest_cancel_function *cancel_fn,
                         MPIX_Grequest_poll_function *poll_fn,
                         MPIX_Grequest_wait_function *wait_fn, void *extra_state,
                         MPI_Request *request) ;


























































































void dlacpy_(char*, int*, int*, double*, int*, double*, int*);
void dgemm_(char*, char*, int*, int*, int*, double*, double*, int*, double*, int*, double*, double*, int*); 


void slacpy_(char*, int*, int*, float*, int*, float*, int*);
void sgemm_(char*, char*, int*, int*, int*, float*, float*, int*, float*, int*, float*, float*, int*); 




void zlacpy_(char*, int*, int*, double _Complex*, int*, double _Complex*, int*);
void zgemm_(char*, char*, int*, int*, int*, double _Complex*, double _Complex*, int*, double _Complex*, int*, double _Complex*, double _Complex*, int*); 


void clacpy_(char*, int*, int*, float _Complex*, int*, float _Complex*, int*);
void cgemm_(char*, char*, int*, int*, int*, float _Complex*, float _Complex*, int*, float _Complex*, int*, float _Complex*, float _Complex*, int*); 




int numroc_(int*, int*, int*, int*, int*);


void pdlacpy_(char*, int*, int*, double*, int*, int*, int*, double*, int*, int*, int*);
void pdtran_(int*, int*, double*, double*, int*, int*, int*, double*, double*, int*, int*, int*);


void pslacpy_(char*, int*, int*, float*, int*, int*, int*, float*, int*, int*, int*);
void pstran_(int*, int*, float*, float*, int*, int*, int*, float*, float*, int*, int*, int*);



void pzlacpy_(char*, int*, int*, double _Complex*, int*, int*, int*, double _Complex*, int*, int*, int*);
void pztranc_(int*, int*, double _Complex*, double _Complex*, int*, int*, int*, double _Complex*, double _Complex*, int*, int*, int*);


void pclacpy_(char*, int*, int*, float _Complex*, int*, int*, int*, float _Complex*, int*, int*, int*);
void pctranc_(int*, int*, float _Complex*, float _Complex*, int*, int*, int*, float _Complex*, float _Complex*, int*, int*, int*);



void cannons_reduction_d(double* A, double* U, int np_rows, int np_cols, int my_prow, int my_pcol,
                         int* a_desc, double *Res, int ToStore, MPI_Comm row_comm, MPI_Comm col_comm)
{
   
      
      
   
      
   
   
  
   int na, nblk, i, j, Size_send_A, Size_receive_A, Size_send_U, Size_receive_U, Buf_rows, Buf_cols, where_to_send_A, from_where_to_receive_A, where_to_send_U, from_where_to_receive_U, last_proc_row, last_proc_col, cols_in_buffer_A, rows_in_buffer_A, intNumber;
   double *Buf_to_send_A, *Buf_to_receive_A, *Buf_to_send_U, *Buf_to_receive_U, *data_ptr, *Buf_A, *Buf_pos, *U_local_start, *Res_ptr, *M, *M_T, *A_local_start, *U_local_start_curr, *U_stored, *CopyTo, *CopyFrom, *U_to_calc;
   int ratio, num_of_iters, cols_in_buffer, rows_in_block, rows_in_buffer, curr_col_loc, cols_in_block, curr_col_glob, curr_row_loc, Size_receive_A_now, Nb, owner, cols_in_buffer_A_now;
   int Size_receive_A_nowMPI, Size_receive_AMPI, Size_receive_UMPI;

   int  row_of_origin_U, rows_in_block_U, num_of_blocks_in_U_buffer, k, startPos, cols_in_buffer_U, rows_in_buffer_U, col_of_origin_A, curr_row_loc_res, curr_row_loc_A, curr_col_glob_res; 
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
   na_rows = numroc_(&na, &nblk, &my_prow, &zero, &np_rows);
   na_cols = numroc_(&na, &nblk, &my_pcol, &zero, &np_cols); 
   
   
   
   
   
   
   

   if (np_cols%np_rows != 0)
   {
      
      
      return;
   }
   if (np_cols < np_rows != 0)
   {
      
      
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
         Buf_rows = na_rows + 1;   
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
      dlacpy_("A", &na_rows, &na_cols, A, &na_rows, Buf_to_send_A, &na_rows);   
   Size_receive_A = 0; 
   
   
   for(i = 0; i < ratio; i++)
   {
      where_to_send_A = (my_pcol - my_prow - i*np_rows + np_cols)%np_cols;                
      from_where_to_receive_A = (my_pcol + my_prow + i*np_rows)%np_cols;
      
      
      if(ratio != 1)   
      {
         if(where_to_send_A != my_pcol)
         {
           MPI_Sendrecv(Buf_to_send_A, (int) (na_cols*na_rows), ((MPI_Datatype)0x4c00080b),(int) where_to_send_A, (int) zero, Buf_A, (int) (na_rows*Buf_cols), ((MPI_Datatype)0x4c00080b), (int) from_where_to_receive_A, (int) zero, row_comm, &status);
           MPI_Get_count(&status, ((MPI_Datatype)0x4c00080b), &Size_receive_A_nowMPI);
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
            dlacpy_("A", &na_rows, &width, CopyFrom, &na_rows, CopyTo, &na_rows);
            CopyTo = CopyTo + na_rows*nblk*ratio; 
            CopyFrom = CopyFrom + na_rows*nblk; 
         }
      }
      else  
         if(my_prow > 0)
         {
            dlacpy_("A", &na_rows, &na_cols, A, &na_rows, Buf_to_send_A, &na_rows);   
            MPI_Sendrecv(Buf_to_send_A, (int) (na_cols*na_rows), ((MPI_Datatype)0x4c00080b), (int) where_to_send_A, (int) zero, Buf_to_receive_A, (int) (na_rows*Buf_cols), ((MPI_Datatype)0x4c00080b), (int) from_where_to_receive_A, (int) zero, row_comm, &status);
            MPI_Get_count(&status, ((MPI_Datatype)0x4c00080b), &Size_receive_AMPI);
            Size_receive_A = (int) Size_receive_AMPI;
	    Size_receive_A = Size_receive_A/na_rows;       
         }
         else
         {
            dlacpy_("A", &na_rows, &na_cols, A, &na_rows, Buf_to_receive_A, &na_rows);   
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
         dlacpy_("A", &rows_in_block, &cols_in_block, data_ptr, &na_rows, Buf_pos, &rows_in_block);     
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
      
      MPI_Sendrecv(Buf_to_send_U, (int) Size_send_U, ((MPI_Datatype)0x4c00080b), (int) where_to_send_U, (int) zero, Buf_to_receive_U, (int) (Buf_rows*na_cols), ((MPI_Datatype)0x4c00080b), (int) from_where_to_receive_U, (int) zero, col_comm, &status); 
      MPI_Get_count(&status, ((MPI_Datatype)0x4c00080b), &Size_receive_UMPI); 
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
      MPI_Isend(Buf_to_send_A, (int) (Size_send_A*na_rows), ((MPI_Datatype)0x4c00080b), (int) where_to_send_A, (int) zero, row_comm, &request_A_Send); 
      MPI_Irecv(Buf_to_receive_A, (int) (Buf_cols*na_rows*ratio), ((MPI_Datatype)0x4c00080b), (int) from_where_to_receive_A, (int) zero, row_comm, &request_A_Recv);
         
      
      Size_send_U = Size_receive_U; 
      MPI_Isend(Buf_to_send_U, (int) Size_send_U, ((MPI_Datatype)0x4c00080b), (int) where_to_send_U, (int) zero, col_comm, &request_U_Send); 
      MPI_Irecv(Buf_to_receive_U, (int) (Buf_rows*na_cols), ((MPI_Datatype)0x4c00080b), (int) from_where_to_receive_U, (int) zero, col_comm, &request_U_Recv); 
      
      
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

      MPI_Get_count(&status, ((MPI_Datatype)0x4c00080b), &Size_receive_AMPI); 
      Size_receive_A = (int) Size_receive_AMPI;
      Size_receive_A = Size_receive_A / na_rows;
      
      
      MPI_Wait(&request_U_Send, &status);
      MPI_Wait(&request_U_Recv, &status);
      MPI_Get_count(&status, ((MPI_Datatype)0x4c00080b), &Size_receive_UMPI); 
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
   
   
   
   pdtran_(&na, &na, &done, M, &one, &one, a_desc, &dzero, M_T, &one, &one, a_desc);     
 
   
           
   
   
   
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
         dlacpy_("A", &rows_in_block, &cols_in_block, A_local_start, &na_rows, Buf_pos, &rows_in_block);     
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
            MPI_Sendrecv(Buf_to_send_A, (int) Size_send_A, ((MPI_Datatype)0x4c00080b), (int) where_to_send_A, (int) zero, Buf_A, (int) Size_U_stored, ((MPI_Datatype)0x4c00080b), (int) from_where_to_receive_A, (int) zero, row_comm, &status);
            MPI_Get_count(&status, ((MPI_Datatype)0x4c00080b), &Size_receive_A_nowMPI);
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
               
            dlacpy_("A", &rows_in_block, &cols_in_block, CopyFrom, &rows_in_block, CopyTo, &rows_in_block);

            CopyFrom = CopyFrom + rows_in_block*cols_in_block; 
            CopyTo = CopyTo + nblk*(ratio*rows_in_block - nblk*(ratio-1)*ratio/2);  
            rows_in_block = rows_in_block - ratio*nblk;     
         }
      }
      else    
      {
         if(my_prow > 0)
         {
            MPI_Sendrecv(Buf_to_send_A, (int) Size_send_A, ((MPI_Datatype)0x4c00080b), (int) where_to_send_A, (int) zero, Buf_to_receive_A, (int) Size_U_stored, ((MPI_Datatype)0x4c00080b), (int) from_where_to_receive_A, (int) zero, row_comm, &status);
            MPI_Get_count(&status, ((MPI_Datatype)0x4c00080b), &Size_receive_AMPI);
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
      MPI_Isend(Buf_to_send_A, (int) Size_send_A, ((MPI_Datatype)0x4c00080b), (int) where_to_send_A, (int) zero, row_comm, &request_A_Send); 
      MPI_Irecv(Buf_to_receive_A, (int) (ratio*Size_U_stored), ((MPI_Datatype)0x4c00080b), (int) from_where_to_receive_A, (int) zero, row_comm, &request_A_Recv);
         
      
      Size_send_U = Size_receive_U; 
      if (j > ToStore)
      {
         if(j > ToStore + 1)
         {
            MPI_Isend(Buf_to_send_U, (int) Size_send_U, ((MPI_Datatype)0x4c00080b), (int) where_to_send_U, (int) zero, col_comm, &request_U_Send);
            U_to_calc = Buf_to_send_U;
         }
         else {
	    MPI_Isend(U_to_calc, (int) Size_send_U, ((MPI_Datatype)0x4c00080b), (int) where_to_send_U, (int) zero, col_comm, &request_U_Send);
	 }
         MPI_Irecv(Buf_to_receive_U, (int) Size_U_stored, ((MPI_Datatype)0x4c00080b), (int) from_where_to_receive_U, (int) zero, col_comm, &request_U_Recv);	 
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
      MPI_Get_count(&status, ((MPI_Datatype)0x4c00080b), &Size_receive_AMPI); 
      Size_receive_A = (int) Size_receive_AMPI;
      
      if (j <= ToStore)
      {
         U_to_calc = &U_stored[Curr_pos_in_U_stored];
         Curr_pos_in_U_stored = Curr_pos_in_U_stored + SizesU[j-1]; 
         Size_receive_U =  SizesU[j-1];
      }
      else
      {
         MPI_Wait(&request_U_Send, &status);
         MPI_Wait(&request_U_Recv, &status);
	 MPI_Get_count(&status, ((MPI_Datatype)0x4c00080b), &Size_receive_UMPI); 
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
   
   pdtran_(&na, &na, &done, Res, &one, &one, a_desc, &dzero, M, &one, &one, a_desc);
   pdlacpy_("U", &na, &na, M, &one, &one, a_desc, Res, &one, &one, a_desc);
      

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

  MPI_Comm c_row_comm = (MPI_Comm)(row_comm);
  MPI_Comm c_col_comm = (MPI_Comm)(col_comm);


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
   
   na_rows = numroc_(&na, &nblk, &my_prow, &zero, &np_rows);
   na_cols = numroc_(&na, &nblk, &my_pcol, &zero, &np_cols);
   nb_cols = numroc_(&nb, &nblk, &my_pcol, &zero, &np_cols);
   
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
         dlacpy_("A", &rows_in_block, &cols_in_block, double_ptr, &na_rows, Buf_pos, &rows_in_block);     
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
            MPI_Sendrecv(Buf_to_send_U, (int) Size_send_U, ((MPI_Datatype)0x4c00080b), (int) where_to_send_U, 0, Buf_U, (int) Size_U_stored, ((MPI_Datatype)0x4c00080b), (int) from_where_to_receive_U, 0, row_comm, &status);
            MPI_Get_count(&status, ((MPI_Datatype)0x4c00080b), &Size_receive_U_nowMPI);
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
               
            dlacpy_("A", &rows_in_block, &cols_in_block, CopyFrom, &rows_in_block, CopyTo, &rows_in_block);

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
            MPI_Sendrecv(Buf_to_send_U, (int) Size_send_U, ((MPI_Datatype)0x4c00080b), (int) where_to_send_U, 0, Buf_to_receive_U, (int) Size_U_stored, ((MPI_Datatype)0x4c00080b), (int) from_where_to_receive_U, 0, row_comm, &status);
            MPI_Get_count(&status, ((MPI_Datatype)0x4c00080b), &Size_receive_UMPI);
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
         
         dlacpy_("A", &na_rows, &nb_cols, B, &na_rows, Buf_to_send_B, &na_rows);
         MPI_Sendrecv(Buf_to_send_B, (int) (nb_cols*na_rows), ((MPI_Datatype)0x4c00080b), (int) where_to_send_B, 0, Buf_to_receive_B, (int) (nb_cols*Buf_rows), ((MPI_Datatype)0x4c00080b), (int) from_where_to_receive_B, 0, col_comm, &status); 
         MPI_Get_count(&status, ((MPI_Datatype)0x4c00080b), &Size_receive_BMPI); 
         Size_receive_B = (int) Size_receive_BMPI;
         Size_receive_B = Size_receive_B/nb_cols;    
	 
      }
      else
      {
         dlacpy_("A", &na_rows, &nb_cols, B, &na_rows, Buf_to_receive_B, &na_rows); 
         Size_receive_B = na_rows;
      }
   }
   else
   {
      dlacpy_("A", &na_rows, &nb_cols, B, &na_rows, Buf_to_receive_B, &na_rows);        
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
        
      
      MPI_Isend(Buf_to_send_U, (int) Size_send_U, ((MPI_Datatype)0x4c00080b), (int) where_to_send_U, 0, row_comm, &request_U_Send); 
      MPI_Irecv(Buf_to_receive_U, (int) (ratio*Size_U_stored), ((MPI_Datatype)0x4c00080b), (int) from_where_to_receive_U, 0, row_comm, &request_U_Recv);      
      
      MPI_Isend(Buf_to_send_B, (int) (Size_send_B*nb_cols), ((MPI_Datatype)0x4c00080b), (int) where_to_send_B, 0, col_comm, &request_B_Send); 
      MPI_Irecv(Buf_to_receive_B, (int) (Buf_rows*nb_cols), ((MPI_Datatype)0x4c00080b), (int) from_where_to_receive_B, 0, col_comm, &request_B_Recv);      
      
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
      MPI_Get_count(&status, ((MPI_Datatype)0x4c00080b), &Size_receive_UMPI); 
      Size_receive_U = (int) Size_receive_UMPI;

      MPI_Wait(&request_B_Send, &status);
      MPI_Wait(&request_B_Recv, &status);
      MPI_Get_count(&status, ((MPI_Datatype)0x4c00080b), &Size_receive_BMPI); 
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

  MPI_Comm c_row_comm = (MPI_Comm)(row_comm);
  MPI_Comm c_col_comm = (MPI_Comm)(col_comm);

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























































































void dlacpy_(char*, int*, int*, double*, int*, double*, int*);
void dgemm_(char*, char*, int*, int*, int*, double*, double*, int*, double*, int*, double*, double*, int*); 


void slacpy_(char*, int*, int*, float*, int*, float*, int*);
void sgemm_(char*, char*, int*, int*, int*, float*, float*, int*, float*, int*, float*, float*, int*); 




void zlacpy_(char*, int*, int*, double _Complex*, int*, double _Complex*, int*);
void zgemm_(char*, char*, int*, int*, int*, double _Complex*, double _Complex*, int*, double _Complex*, int*, double _Complex*, double _Complex*, int*); 


void clacpy_(char*, int*, int*, float _Complex*, int*, float _Complex*, int*);
void cgemm_(char*, char*, int*, int*, int*, float _Complex*, float _Complex*, int*, float _Complex*, int*, float _Complex*, float _Complex*, int*); 




int numroc_(int*, int*, int*, int*, int*);


void pdlacpy_(char*, int*, int*, double*, int*, int*, int*, double*, int*, int*, int*);
void pdtran_(int*, int*, double*, double*, int*, int*, int*, double*, double*, int*, int*, int*);


void pslacpy_(char*, int*, int*, float*, int*, int*, int*, float*, int*, int*, int*);
void pstran_(int*, int*, float*, float*, int*, int*, int*, float*, float*, int*, int*, int*);



void pzlacpy_(char*, int*, int*, double _Complex*, int*, int*, int*, double _Complex*, int*, int*, int*);
void pztranc_(int*, int*, double _Complex*, double _Complex*, int*, int*, int*, double _Complex*, double _Complex*, int*, int*, int*);


void pclacpy_(char*, int*, int*, float _Complex*, int*, int*, int*, float _Complex*, int*, int*, int*);
void pctranc_(int*, int*, float _Complex*, float _Complex*, int*, int*, int*, float _Complex*, float _Complex*, int*, int*, int*);



void cannons_reduction_f(float* A, float* U, int np_rows, int np_cols, int my_prow, int my_pcol,
                         int* a_desc, float *Res, int ToStore, MPI_Comm row_comm, MPI_Comm col_comm)
{
   
      
      
   
      
   
   
  
   int na, nblk, i, j, Size_send_A, Size_receive_A, Size_send_U, Size_receive_U, Buf_rows, Buf_cols, where_to_send_A, from_where_to_receive_A, where_to_send_U, from_where_to_receive_U, last_proc_row, last_proc_col, cols_in_buffer_A, rows_in_buffer_A, intNumber;
   float *Buf_to_send_A, *Buf_to_receive_A, *Buf_to_send_U, *Buf_to_receive_U, *data_ptr, *Buf_A, *Buf_pos, *U_local_start, *Res_ptr, *M, *M_T, *A_local_start, *U_local_start_curr, *U_stored, *CopyTo, *CopyFrom, *U_to_calc;
   int ratio, num_of_iters, cols_in_buffer, rows_in_block, rows_in_buffer, curr_col_loc, cols_in_block, curr_col_glob, curr_row_loc, Size_receive_A_now, Nb, owner, cols_in_buffer_A_now;
   int Size_receive_A_nowMPI, Size_receive_AMPI, Size_receive_UMPI;

   int  row_of_origin_U, rows_in_block_U, num_of_blocks_in_U_buffer, k, startPos, cols_in_buffer_U, rows_in_buffer_U, col_of_origin_A, curr_row_loc_res, curr_row_loc_A, curr_col_glob_res; 
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
   na_rows = numroc_(&na, &nblk, &my_prow, &zero, &np_rows);
   na_cols = numroc_(&na, &nblk, &my_pcol, &zero, &np_cols); 
   
   
   
   
   
   
   

   if (np_cols%np_rows != 0)
   {
      
      
      return;
   }
   if (np_cols < np_rows != 0)
   {
      
      
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
         Buf_rows = na_rows + 1;   
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
      slacpy_("A", &na_rows, &na_cols, A, &na_rows, Buf_to_send_A, &na_rows);   
   Size_receive_A = 0; 
   
   
   for(i = 0; i < ratio; i++)
   {
      where_to_send_A = (my_pcol - my_prow - i*np_rows + np_cols)%np_cols;                
      from_where_to_receive_A = (my_pcol + my_prow + i*np_rows)%np_cols;
      
      
      if(ratio != 1)   
      {
         if(where_to_send_A != my_pcol)
         {
           MPI_Sendrecv(Buf_to_send_A, (int) (na_cols*na_rows), ((MPI_Datatype)0x4c00040a),(int) where_to_send_A, (int) zero, Buf_A, (int) (na_rows*Buf_cols), ((MPI_Datatype)0x4c00040a), (int) from_where_to_receive_A, (int) zero, row_comm, &status);
           MPI_Get_count(&status, ((MPI_Datatype)0x4c00040a), &Size_receive_A_nowMPI);
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
            slacpy_("A", &na_rows, &width, CopyFrom, &na_rows, CopyTo, &na_rows);
            CopyTo = CopyTo + na_rows*nblk*ratio; 
            CopyFrom = CopyFrom + na_rows*nblk; 
         }
      }
      else  
         if(my_prow > 0)
         {
            slacpy_("A", &na_rows, &na_cols, A, &na_rows, Buf_to_send_A, &na_rows);   
            MPI_Sendrecv(Buf_to_send_A, (int) (na_cols*na_rows), ((MPI_Datatype)0x4c00040a), (int) where_to_send_A, (int) zero, Buf_to_receive_A, (int) (na_rows*Buf_cols), ((MPI_Datatype)0x4c00040a), (int) from_where_to_receive_A, (int) zero, row_comm, &status);
            MPI_Get_count(&status, ((MPI_Datatype)0x4c00040a), &Size_receive_AMPI);
            Size_receive_A = (int) Size_receive_AMPI;
	    Size_receive_A = Size_receive_A/na_rows;       
         }
         else
         {
            slacpy_("A", &na_rows, &na_cols, A, &na_rows, Buf_to_receive_A, &na_rows);   
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
         slacpy_("A", &rows_in_block, &cols_in_block, data_ptr, &na_rows, Buf_pos, &rows_in_block);     
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
      
      MPI_Sendrecv(Buf_to_send_U, (int) Size_send_U, ((MPI_Datatype)0x4c00040a), (int) where_to_send_U, (int) zero, Buf_to_receive_U, (int) (Buf_rows*na_cols), ((MPI_Datatype)0x4c00040a), (int) from_where_to_receive_U, (int) zero, col_comm, &status); 
      MPI_Get_count(&status, ((MPI_Datatype)0x4c00040a), &Size_receive_UMPI); 
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
      MPI_Isend(Buf_to_send_A, (int) (Size_send_A*na_rows), ((MPI_Datatype)0x4c00040a), (int) where_to_send_A, (int) zero, row_comm, &request_A_Send); 
      MPI_Irecv(Buf_to_receive_A, (int) (Buf_cols*na_rows*ratio), ((MPI_Datatype)0x4c00040a), (int) from_where_to_receive_A, (int) zero, row_comm, &request_A_Recv);
         
      
      Size_send_U = Size_receive_U; 
      MPI_Isend(Buf_to_send_U, (int) Size_send_U, ((MPI_Datatype)0x4c00040a), (int) where_to_send_U, (int) zero, col_comm, &request_U_Send); 
      MPI_Irecv(Buf_to_receive_U, (int) (Buf_rows*na_cols), ((MPI_Datatype)0x4c00040a), (int) from_where_to_receive_U, (int) zero, col_comm, &request_U_Recv); 
      
      
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

      MPI_Get_count(&status, ((MPI_Datatype)0x4c00040a), &Size_receive_AMPI); 
      Size_receive_A = (int) Size_receive_AMPI;
      Size_receive_A = Size_receive_A / na_rows;
      
      
      MPI_Wait(&request_U_Send, &status);
      MPI_Wait(&request_U_Recv, &status);
      MPI_Get_count(&status, ((MPI_Datatype)0x4c00040a), &Size_receive_UMPI); 
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
   
   
   
   pstran_(&na, &na, &done, M, &one, &one, a_desc, &dzero, M_T, &one, &one, a_desc);     
 
   
           
   
   
   
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
         slacpy_("A", &rows_in_block, &cols_in_block, A_local_start, &na_rows, Buf_pos, &rows_in_block);     
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
            MPI_Sendrecv(Buf_to_send_A, (int) Size_send_A, ((MPI_Datatype)0x4c00040a), (int) where_to_send_A, (int) zero, Buf_A, (int) Size_U_stored, ((MPI_Datatype)0x4c00040a), (int) from_where_to_receive_A, (int) zero, row_comm, &status);
            MPI_Get_count(&status, ((MPI_Datatype)0x4c00040a), &Size_receive_A_nowMPI);
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
               
            slacpy_("A", &rows_in_block, &cols_in_block, CopyFrom, &rows_in_block, CopyTo, &rows_in_block);

            CopyFrom = CopyFrom + rows_in_block*cols_in_block; 
            CopyTo = CopyTo + nblk*(ratio*rows_in_block - nblk*(ratio-1)*ratio/2);  
            rows_in_block = rows_in_block - ratio*nblk;     
         }
      }
      else    
      {
         if(my_prow > 0)
         {
            MPI_Sendrecv(Buf_to_send_A, (int) Size_send_A, ((MPI_Datatype)0x4c00040a), (int) where_to_send_A, (int) zero, Buf_to_receive_A, (int) Size_U_stored, ((MPI_Datatype)0x4c00040a), (int) from_where_to_receive_A, (int) zero, row_comm, &status);
            MPI_Get_count(&status, ((MPI_Datatype)0x4c00040a), &Size_receive_AMPI);
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
      MPI_Isend(Buf_to_send_A, (int) Size_send_A, ((MPI_Datatype)0x4c00040a), (int) where_to_send_A, (int) zero, row_comm, &request_A_Send); 
      MPI_Irecv(Buf_to_receive_A, (int) (ratio*Size_U_stored), ((MPI_Datatype)0x4c00040a), (int) from_where_to_receive_A, (int) zero, row_comm, &request_A_Recv);
         
      
      Size_send_U = Size_receive_U; 
      if (j > ToStore)
      {
         if(j > ToStore + 1)
         {
            MPI_Isend(Buf_to_send_U, (int) Size_send_U, ((MPI_Datatype)0x4c00040a), (int) where_to_send_U, (int) zero, col_comm, &request_U_Send);
            U_to_calc = Buf_to_send_U;
         }
         else {
	    MPI_Isend(U_to_calc, (int) Size_send_U, ((MPI_Datatype)0x4c00040a), (int) where_to_send_U, (int) zero, col_comm, &request_U_Send);
	 }
         MPI_Irecv(Buf_to_receive_U, (int) Size_U_stored, ((MPI_Datatype)0x4c00040a), (int) from_where_to_receive_U, (int) zero, col_comm, &request_U_Recv);	 
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
      MPI_Get_count(&status, ((MPI_Datatype)0x4c00040a), &Size_receive_AMPI); 
      Size_receive_A = (int) Size_receive_AMPI;
      
      if (j <= ToStore)
      {
         U_to_calc = &U_stored[Curr_pos_in_U_stored];
         Curr_pos_in_U_stored = Curr_pos_in_U_stored + SizesU[j-1]; 
         Size_receive_U =  SizesU[j-1];
      }
      else
      {
         MPI_Wait(&request_U_Send, &status);
         MPI_Wait(&request_U_Recv, &status);
	 MPI_Get_count(&status, ((MPI_Datatype)0x4c00040a), &Size_receive_UMPI); 
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
   
   pstran_(&na, &na, &done, Res, &one, &one, a_desc, &dzero, M, &one, &one, a_desc);
   pslacpy_("U", &na, &na, M, &one, &one, a_desc, Res, &one, &one, a_desc);
      

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

  MPI_Comm c_row_comm = (MPI_Comm)(row_comm);
  MPI_Comm c_col_comm = (MPI_Comm)(col_comm);


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
   
   na_rows = numroc_(&na, &nblk, &my_prow, &zero, &np_rows);
   na_cols = numroc_(&na, &nblk, &my_pcol, &zero, &np_cols);
   nb_cols = numroc_(&nb, &nblk, &my_pcol, &zero, &np_cols);
   
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
         slacpy_("A", &rows_in_block, &cols_in_block, double_ptr, &na_rows, Buf_pos, &rows_in_block);     
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
            MPI_Sendrecv(Buf_to_send_U, (int) Size_send_U, ((MPI_Datatype)0x4c00040a), (int) where_to_send_U, 0, Buf_U, (int) Size_U_stored, ((MPI_Datatype)0x4c00040a), (int) from_where_to_receive_U, 0, row_comm, &status);
            MPI_Get_count(&status, ((MPI_Datatype)0x4c00040a), &Size_receive_U_nowMPI);
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
               
            slacpy_("A", &rows_in_block, &cols_in_block, CopyFrom, &rows_in_block, CopyTo, &rows_in_block);

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
            MPI_Sendrecv(Buf_to_send_U, (int) Size_send_U, ((MPI_Datatype)0x4c00040a), (int) where_to_send_U, 0, Buf_to_receive_U, (int) Size_U_stored, ((MPI_Datatype)0x4c00040a), (int) from_where_to_receive_U, 0, row_comm, &status);
            MPI_Get_count(&status, ((MPI_Datatype)0x4c00040a), &Size_receive_UMPI);
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
         
         slacpy_("A", &na_rows, &nb_cols, B, &na_rows, Buf_to_send_B, &na_rows);
         MPI_Sendrecv(Buf_to_send_B, (int) (nb_cols*na_rows), ((MPI_Datatype)0x4c00040a), (int) where_to_send_B, 0, Buf_to_receive_B, (int) (nb_cols*Buf_rows), ((MPI_Datatype)0x4c00040a), (int) from_where_to_receive_B, 0, col_comm, &status); 
         MPI_Get_count(&status, ((MPI_Datatype)0x4c00040a), &Size_receive_BMPI); 
         Size_receive_B = (int) Size_receive_BMPI;
         Size_receive_B = Size_receive_B/nb_cols;    
	 
      }
      else
      {
         slacpy_("A", &na_rows, &nb_cols, B, &na_rows, Buf_to_receive_B, &na_rows); 
         Size_receive_B = na_rows;
      }
   }
   else
   {
      slacpy_("A", &na_rows, &nb_cols, B, &na_rows, Buf_to_receive_B, &na_rows);        
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
        
      
      MPI_Isend(Buf_to_send_U, (int) Size_send_U, ((MPI_Datatype)0x4c00040a), (int) where_to_send_U, 0, row_comm, &request_U_Send); 
      MPI_Irecv(Buf_to_receive_U, (int) (ratio*Size_U_stored), ((MPI_Datatype)0x4c00040a), (int) from_where_to_receive_U, 0, row_comm, &request_U_Recv);      
      
      MPI_Isend(Buf_to_send_B, (int) (Size_send_B*nb_cols), ((MPI_Datatype)0x4c00040a), (int) where_to_send_B, 0, col_comm, &request_B_Send); 
      MPI_Irecv(Buf_to_receive_B, (int) (Buf_rows*nb_cols), ((MPI_Datatype)0x4c00040a), (int) from_where_to_receive_B, 0, col_comm, &request_B_Recv);      
      
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
      MPI_Get_count(&status, ((MPI_Datatype)0x4c00040a), &Size_receive_UMPI); 
      Size_receive_U = (int) Size_receive_UMPI;

      MPI_Wait(&request_B_Send, &status);
      MPI_Wait(&request_B_Recv, &status);
      MPI_Get_count(&status, ((MPI_Datatype)0x4c00040a), &Size_receive_BMPI); 
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

  MPI_Comm c_row_comm = (MPI_Comm)(row_comm);
  MPI_Comm c_col_comm = (MPI_Comm)(col_comm);

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






















































































void dlacpy_(char*, int*, int*, double*, int*, double*, int*);
void dgemm_(char*, char*, int*, int*, int*, double*, double*, int*, double*, int*, double*, double*, int*); 


void slacpy_(char*, int*, int*, float*, int*, float*, int*);
void sgemm_(char*, char*, int*, int*, int*, float*, float*, int*, float*, int*, float*, float*, int*); 




void zlacpy_(char*, int*, int*, double _Complex*, int*, double _Complex*, int*);
void zgemm_(char*, char*, int*, int*, int*, double _Complex*, double _Complex*, int*, double _Complex*, int*, double _Complex*, double _Complex*, int*); 


void clacpy_(char*, int*, int*, float _Complex*, int*, float _Complex*, int*);
void cgemm_(char*, char*, int*, int*, int*, float _Complex*, float _Complex*, int*, float _Complex*, int*, float _Complex*, float _Complex*, int*); 




int numroc_(int*, int*, int*, int*, int*);


void pdlacpy_(char*, int*, int*, double*, int*, int*, int*, double*, int*, int*, int*);
void pdtran_(int*, int*, double*, double*, int*, int*, int*, double*, double*, int*, int*, int*);


void pslacpy_(char*, int*, int*, float*, int*, int*, int*, float*, int*, int*, int*);
void pstran_(int*, int*, float*, float*, int*, int*, int*, float*, float*, int*, int*, int*);



void pzlacpy_(char*, int*, int*, double _Complex*, int*, int*, int*, double _Complex*, int*, int*, int*);
void pztranc_(int*, int*, double _Complex*, double _Complex*, int*, int*, int*, double _Complex*, double _Complex*, int*, int*, int*);


void pclacpy_(char*, int*, int*, float _Complex*, int*, int*, int*, float _Complex*, int*, int*, int*);
void pctranc_(int*, int*, float _Complex*, float _Complex*, int*, int*, int*, float _Complex*, float _Complex*, int*, int*, int*);



void cannons_reduction_dc(double _Complex* A, double _Complex* U, int np_rows, int np_cols, int my_prow, int my_pcol,
                         int* a_desc, double _Complex *Res, int ToStore, MPI_Comm row_comm, MPI_Comm col_comm)
{
   
      
      
   
      
   
   
  
   int na, nblk, i, j, Size_send_A, Size_receive_A, Size_send_U, Size_receive_U, Buf_rows, Buf_cols, where_to_send_A, from_where_to_receive_A, where_to_send_U, from_where_to_receive_U, last_proc_row, last_proc_col, cols_in_buffer_A, rows_in_buffer_A, intNumber;
   double _Complex *Buf_to_send_A, *Buf_to_receive_A, *Buf_to_send_U, *Buf_to_receive_U, *data_ptr, *Buf_A, *Buf_pos, *U_local_start, *Res_ptr, *M, *M_T, *A_local_start, *U_local_start_curr, *U_stored, *CopyTo, *CopyFrom, *U_to_calc;
   int ratio, num_of_iters, cols_in_buffer, rows_in_block, rows_in_buffer, curr_col_loc, cols_in_block, curr_col_glob, curr_row_loc, Size_receive_A_now, Nb, owner, cols_in_buffer_A_now;
   int Size_receive_A_nowMPI, Size_receive_AMPI, Size_receive_UMPI;

   int  row_of_origin_U, rows_in_block_U, num_of_blocks_in_U_buffer, k, startPos, cols_in_buffer_U, rows_in_buffer_U, col_of_origin_A, curr_row_loc_res, curr_row_loc_A, curr_col_glob_res; 
   int curr_col_loc_res, curr_col_loc_buf, proc_row_curr, curr_col_loc_U, A_local_index, LDA_A, LDA_A_new, index_row_A_for_LDA, ii, rows_in_block_U_curr, width, row_origin_U, rows_in_block_A, cols_in_buffer_A_my_initial, rows_in_buffer_A_my_initial, proc_col_min;
   int *SizesU;
   int Size_U_skewed, Size_U_stored, Curr_pos_in_U_stored, rows_in_buffer_A_now;
   double _Complex done = 1.0;
   double _Complex dzero = 0.0;
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
   na_rows = numroc_(&na, &nblk, &my_prow, &zero, &np_rows);
   na_cols = numroc_(&na, &nblk, &my_pcol, &zero, &np_cols); 
   
   
   
   
   
   
   

   if (np_cols%np_rows != 0)
   {
      
      
      return;
   }
   if (np_cols < np_rows != 0)
   {
      
      
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
         Buf_rows = na_rows + 1;   
      else
         Buf_rows = na_rows + nblk;      
   else
      if(my_prow < last_proc_row)
         Buf_rows = na_rows;
      else if(my_prow > last_proc_row)
         Buf_rows = na_rows + nblk; 
      else  
         Buf_rows = na_rows + nblk - na_rows%nblk;  
      
   intNumber = ceil((double _Complex)na/(double _Complex)(np_cols*nblk));   
   Size_U_stored = ratio*nblk*nblk*intNumber*(intNumber+1)/2 + 2;   
   
   U_stored = malloc((Size_U_stored*(ToStore+1))*sizeof(double _Complex));
   SizesU = malloc(ToStore*sizeof(int));  
   Buf_to_send_A = malloc(ratio*Buf_cols*Buf_rows*sizeof(double _Complex));
   Buf_to_receive_A = malloc(ratio*Buf_cols*Buf_rows*sizeof(double _Complex));
   Buf_to_send_U = malloc(Size_U_stored*sizeof(double _Complex));
   Buf_to_receive_U = malloc(Size_U_stored*sizeof(double _Complex));
   if(ratio != 1)
      Buf_A = malloc(Buf_cols*Buf_rows*sizeof(double _Complex));   
   M = malloc(na_rows*na_cols*sizeof(double _Complex));
   M_T = malloc(na_rows*na_cols*sizeof(double _Complex));
   for(i = 0; i < na_rows*na_cols; i++)
      M[i] = 0; 
        
   
   
   
   if(ratio != 1)
      zlacpy_("A", &na_rows, &na_cols, A, &na_rows, Buf_to_send_A, &na_rows);   
   Size_receive_A = 0; 
   
   
   for(i = 0; i < ratio; i++)
   {
      where_to_send_A = (my_pcol - my_prow - i*np_rows + np_cols)%np_cols;                
      from_where_to_receive_A = (my_pcol + my_prow + i*np_rows)%np_cols;
      
      
      if(ratio != 1)   
      {
         if(where_to_send_A != my_pcol)
         {
           MPI_Sendrecv(Buf_to_send_A, (int) (na_cols*na_rows), ((MPI_Datatype)1275072546),(int) where_to_send_A, (int) zero, Buf_A, (int) (na_rows*Buf_cols), ((MPI_Datatype)1275072546), (int) from_where_to_receive_A, (int) zero, row_comm, &status);
           MPI_Get_count(&status, ((MPI_Datatype)1275072546), &Size_receive_A_nowMPI);
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
         
         intNumber = ceil((double _Complex)Size_receive_A_now/(double _Complex)nblk);   
         for(j = 0; j < intNumber; j++)
         {
            width = nblk; 
            if(nblk*(j+1) > Size_receive_A_now)
               width = Size_receive_A_now - nblk*j; 
            zlacpy_("A", &na_rows, &width, CopyFrom, &na_rows, CopyTo, &na_rows);
            CopyTo = CopyTo + na_rows*nblk*ratio; 
            CopyFrom = CopyFrom + na_rows*nblk; 
         }
      }
      else  
         if(my_prow > 0)
         {
            zlacpy_("A", &na_rows, &na_cols, A, &na_rows, Buf_to_send_A, &na_rows);   
            MPI_Sendrecv(Buf_to_send_A, (int) (na_cols*na_rows), ((MPI_Datatype)1275072546), (int) where_to_send_A, (int) zero, Buf_to_receive_A, (int) (na_rows*Buf_cols), ((MPI_Datatype)1275072546), (int) from_where_to_receive_A, (int) zero, row_comm, &status);
            MPI_Get_count(&status, ((MPI_Datatype)1275072546), &Size_receive_AMPI);
            Size_receive_A = (int) Size_receive_AMPI;
	    Size_receive_A = Size_receive_A/na_rows;       
         }
         else
         {
            zlacpy_("A", &na_rows, &na_cols, A, &na_rows, Buf_to_receive_A, &na_rows);   
            Size_receive_A = na_cols; 
         }
   }
   
   
     
   
   num_of_iters = ceil((double _Complex)na_cols/(double _Complex)nblk);             
   
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
      rows_in_block = ceil(((double _Complex)(my_pcol + 1) - (double _Complex)my_prow)/(double _Complex)np_rows)*nblk;
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
         zlacpy_("A", &rows_in_block, &cols_in_block, data_ptr, &na_rows, Buf_pos, &rows_in_block);     
         Buf_pos = Buf_pos + rows_in_block*cols_in_block;                         
         Size_send_U = Size_send_U + rows_in_block*cols_in_block; 
      }
      curr_col_loc = curr_col_loc + nblk;      
      rows_in_block = rows_in_block + ratio*nblk;
   }
   rows_in_buffer = rows_in_block - ratio*nblk;    
   *Buf_pos = (double _Complex)rows_in_buffer; 
   Size_send_U = Size_send_U + 1;
   
   
   if(where_to_send_U != my_prow)
   {   
      
      MPI_Sendrecv(Buf_to_send_U, (int) Size_send_U, ((MPI_Datatype)1275072546), (int) where_to_send_U, (int) zero, Buf_to_receive_U, (int) (Buf_rows*na_cols), ((MPI_Datatype)1275072546), (int) from_where_to_receive_U, (int) zero, col_comm, &status); 
      MPI_Get_count(&status, ((MPI_Datatype)1275072546), &Size_receive_UMPI); 
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
      MPI_Isend(Buf_to_send_A, (int) (Size_send_A*na_rows), ((MPI_Datatype)1275072546), (int) where_to_send_A, (int) zero, row_comm, &request_A_Send); 
      MPI_Irecv(Buf_to_receive_A, (int) (Buf_cols*na_rows*ratio), ((MPI_Datatype)1275072546), (int) from_where_to_receive_A, (int) zero, row_comm, &request_A_Recv);
         
      
      Size_send_U = Size_receive_U; 
      MPI_Isend(Buf_to_send_U, (int) Size_send_U, ((MPI_Datatype)1275072546), (int) where_to_send_U, (int) zero, col_comm, &request_U_Send); 
      MPI_Irecv(Buf_to_receive_U, (int) (Buf_rows*na_cols), ((MPI_Datatype)1275072546), (int) from_where_to_receive_U, (int) zero, col_comm, &request_U_Recv); 
      
      
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
    
      num_of_blocks_in_U_buffer = ceil(((double _Complex)cols_in_buffer - (double _Complex)curr_col_loc_buf)/(double _Complex)nblk); 
      
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

      MPI_Get_count(&status, ((MPI_Datatype)1275072546), &Size_receive_AMPI); 
      Size_receive_A = (int) Size_receive_AMPI;
      Size_receive_A = Size_receive_A / na_rows;
      
      
      MPI_Wait(&request_U_Send, &status);
      MPI_Wait(&request_U_Recv, &status);
      MPI_Get_count(&status, ((MPI_Datatype)1275072546), &Size_receive_UMPI); 
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
    
   num_of_blocks_in_U_buffer = ceil(((double _Complex)cols_in_buffer - (double _Complex)curr_col_loc_buf)/(double _Complex)nblk); 
      
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
   
   
   
   pztranc_(&na, &na, &done, M, &one, &one, a_desc, &dzero, M_T, &one, &one, a_desc);     
 
   
           
   
   
   
   if((ratio != 1)||(my_prow != 0))   
      Buf_pos = Buf_to_send_A;     
   else
      Buf_pos = Buf_to_receive_A;  
   
   
   num_of_iters = ceil((double _Complex)na_cols/(double _Complex)nblk);             
   
   cols_in_buffer_A_my_initial = 0;
   Size_send_A = 0; 
   
   if(my_pcol <= my_prow)  
   {
      curr_row_loc = 0;     
      rows_in_buffer_A_my_initial = na_rows;
   }
   else
   {
      curr_row_loc = ceil((double _Complex)(((double _Complex)my_pcol - (double _Complex)my_prow)/(double _Complex)np_rows))*nblk; 
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
         zlacpy_("A", &rows_in_block, &cols_in_block, A_local_start, &na_rows, Buf_pos, &rows_in_block);     
         Buf_pos = Buf_pos + rows_in_block*cols_in_block;
         Size_send_A = Size_send_A + rows_in_block*cols_in_block; 
         cols_in_buffer_A_my_initial = cols_in_buffer_A_my_initial + cols_in_block; 
      }
      curr_row_loc = curr_row_loc + ratio*nblk;
   }
   *Buf_pos = (double _Complex)cols_in_buffer_A_my_initial; 
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
            MPI_Sendrecv(Buf_to_send_A, (int) Size_send_A, ((MPI_Datatype)1275072546), (int) where_to_send_A, (int) zero, Buf_A, (int) Size_U_stored, ((MPI_Datatype)1275072546), (int) from_where_to_receive_A, (int) zero, row_comm, &status);
            MPI_Get_count(&status, ((MPI_Datatype)1275072546), &Size_receive_A_nowMPI);
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
               rows_in_buffer_A_now = na_rows - ceil((double _Complex)(((double _Complex)from_where_to_receive_A - (double _Complex)my_prow)/(double _Complex)np_rows))*nblk; 
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
            
         
         intNumber = ceil((double _Complex)cols_in_buffer_A_now/(double _Complex)nblk);  
         rows_in_block = rows_in_buffer_A_now; 
         for(j = 0; j < intNumber; j++)
         {
            if((j+1)*nblk < cols_in_buffer_A_now)
               cols_in_block = nblk; 
            else
               cols_in_block = cols_in_buffer_A_now - j*nblk;
               
            zlacpy_("A", &rows_in_block, &cols_in_block, CopyFrom, &rows_in_block, CopyTo, &rows_in_block);

            CopyFrom = CopyFrom + rows_in_block*cols_in_block; 
            CopyTo = CopyTo + nblk*(ratio*rows_in_block - nblk*(ratio-1)*ratio/2);  
            rows_in_block = rows_in_block - ratio*nblk;     
         }
      }
      else    
      {
         if(my_prow > 0)
         {
            MPI_Sendrecv(Buf_to_send_A, (int) Size_send_A, ((MPI_Datatype)1275072546), (int) where_to_send_A, (int) zero, Buf_to_receive_A, (int) Size_U_stored, ((MPI_Datatype)1275072546), (int) from_where_to_receive_A, (int) zero, row_comm, &status);
            MPI_Get_count(&status, ((MPI_Datatype)1275072546), &Size_receive_AMPI);
            Size_receive_A = (int) Size_receive_AMPI;

            cols_in_buffer_A = (int)Buf_to_receive_A[Size_receive_A-1];
            if(from_where_to_receive_A <= my_prow)  
            {
               rows_in_buffer_A = na_rows;
            }
            else
            {
               rows_in_buffer_A = na_rows - ceil((double _Complex)(((double _Complex)from_where_to_receive_A - (double _Complex)my_prow)/(double _Complex)np_rows))*nblk; 
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
      MPI_Isend(Buf_to_send_A, (int) Size_send_A, ((MPI_Datatype)1275072546), (int) where_to_send_A, (int) zero, row_comm, &request_A_Send); 
      MPI_Irecv(Buf_to_receive_A, (int) (ratio*Size_U_stored), ((MPI_Datatype)1275072546), (int) from_where_to_receive_A, (int) zero, row_comm, &request_A_Recv);
         
      
      Size_send_U = Size_receive_U; 
      if (j > ToStore)
      {
         if(j > ToStore + 1)
         {
            MPI_Isend(Buf_to_send_U, (int) Size_send_U, ((MPI_Datatype)1275072546), (int) where_to_send_U, (int) zero, col_comm, &request_U_Send);
            U_to_calc = Buf_to_send_U;
         }
         else {
	    MPI_Isend(U_to_calc, (int) Size_send_U, ((MPI_Datatype)1275072546), (int) where_to_send_U, (int) zero, col_comm, &request_U_Send);
	 }
         MPI_Irecv(Buf_to_receive_U, (int) Size_U_stored, ((MPI_Datatype)1275072546), (int) from_where_to_receive_U, (int) zero, col_comm, &request_U_Recv);	 
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
      
      num_of_blocks_in_U_buffer = ceil((double _Complex)((double _Complex)cols_in_buffer_U/(double _Complex)nblk)); 
      if(my_pcol >= row_of_origin_U)    
         rows_in_block_U = ceil(((double _Complex)(my_pcol + 1) - (double _Complex)row_of_origin_U)/(double _Complex)np_rows)*nblk;  
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
 
            
            for (ii = 0; ii < ceil((double _Complex)rows_in_block_U/(double _Complex)nblk); ii++)
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
      MPI_Get_count(&status, ((MPI_Datatype)1275072546), &Size_receive_AMPI); 
      Size_receive_A = (int) Size_receive_AMPI;
      
      if (j <= ToStore)
      {
         U_to_calc = &U_stored[Curr_pos_in_U_stored];
         Curr_pos_in_U_stored = Curr_pos_in_U_stored + SizesU[j-1]; 
         Size_receive_U =  SizesU[j-1];
      }
      else
      {
         MPI_Wait(&request_U_Send, &status);
         MPI_Wait(&request_U_Recv, &status);
	 MPI_Get_count(&status, ((MPI_Datatype)1275072546), &Size_receive_UMPI); 
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
      
   num_of_blocks_in_U_buffer = ceil((double _Complex)((double _Complex)cols_in_buffer_U/(double _Complex)nblk));
   if(my_pcol >= row_of_origin_U)    
      rows_in_block_U = ceil(((double _Complex)(my_pcol + 1) - (double _Complex)row_of_origin_U)/(double _Complex)np_rows)*nblk;  
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

         
         for (ii = 0; ii < ceil((double _Complex)rows_in_block_U/(double _Complex)nblk); ii++)
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
   
   pztranc_(&na, &na, &done, Res, &one, &one, a_desc, &dzero, M, &one, &one, a_desc);
   pzlacpy_("U", &na, &na, M, &one, &one, a_desc, Res, &one, &one, a_desc);
      

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

void cannons_reduction_c_dc(double _Complex* A, double _Complex* U, int local_rowsCast, int local_colsCast,
                         int* a_desc, double _Complex *Res, int ToStore, int row_comm, int col_comm)
{
  int local_rows, local_cols;
  local_rows = (int) local_rowsCast;
  local_cols = (int) local_colsCast;

  MPI_Comm c_row_comm = (MPI_Comm)(row_comm);
  MPI_Comm c_col_comm = (MPI_Comm)(col_comm);


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





























































void cannons_triang_rectangular_dc(double _Complex* U, double _Complex* B, int np_rows, int np_cols, int my_prow, int my_pcol, int* U_desc, int* b_desc, double _Complex *Res, MPI_Comm row_comm, MPI_Comm col_comm)
{
   
   
   
   
   
   
   
   
   
  
   int na, nb, nblk, width, na_rows, na_cols, nb_cols, cols_in_buffer_U_my_initial, cols_in_buffer_U, rows_in_buffer_U, Size_receive_U_now, rows_in_buffer_U_now, cols_in_buffer_U_now, rows_in_buffer_U_my_initial;

   int Size_receive_U_nowMPI, Size_receive_UMPI, Size_receive_BMPI;
   int i, j, Size_send_U, Size_receive_U, Size_send_B, Size_receive_B, intNumber, Buf_rows, Buf_cols_U, Buf_cols_B, curr_rows, num_of_iters, cols_in_buffer, rows_in_block, curr_col_loc, cols_in_block, num_of_blocks_in_U_buffer, col_of_origin_U, b_rows_mult, b_cols_mult; 
   
   double _Complex *Buf_to_send_U, *Buf_to_receive_U, *Buf_to_send_B, *Buf_to_receive_B, *Buf_U, *PosBuff;
  
   int where_to_send_U, from_where_to_receive_U, where_to_send_B, from_where_to_receive_B, last_proc_col_B, last_proc_row_B, n, Size_U_stored, proc_col_min; 
   
   double _Complex *U_local_start, *Buf_pos, *B_local_start, *double_ptr, *CopyTo, *CopyFrom;
   
   int ratio;
   
   MPI_Status status;

   int one = 1;
   int zero = 0; 
   double _Complex done = 1.0;
   double _Complex dzero = 0.0;
      
   na = U_desc[2];
   nblk = U_desc[4]; 
   nb = b_desc[3];
   
   na_rows = numroc_(&na, &nblk, &my_prow, &zero, &np_rows);
   na_cols = numroc_(&na, &nblk, &my_pcol, &zero, &np_cols);
   nb_cols = numroc_(&nb, &nblk, &my_pcol, &zero, &np_cols);
   
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
   
   intNumber = ceil((double _Complex)na/(double _Complex)(np_cols*nblk));   
   Size_U_stored = ratio*nblk*nblk*intNumber*(intNumber+1)/2 + 2;   
   
   Buf_to_send_U = malloc(ratio*Size_U_stored*sizeof(double _Complex));
   Buf_to_receive_U = malloc(ratio*Size_U_stored*sizeof(double _Complex));
   Buf_to_send_B = malloc(Buf_cols_B*Buf_rows*sizeof(double _Complex));
   Buf_to_receive_B = malloc(Buf_cols_B*Buf_rows*sizeof(double _Complex));
   if(ratio != 1)
      Buf_U = malloc(Size_U_stored*sizeof(double _Complex));   
    
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
      
   num_of_iters = ceil((double _Complex)na_cols/(double _Complex)nblk);             
   num_of_iters = num_of_iters - curr_col_loc;   
   curr_col_loc = curr_col_loc*nblk;             

   if(my_pcol >= my_prow )
      rows_in_block = ceil(((double _Complex)(my_pcol + 1) - (double _Complex)my_prow)/(double _Complex)np_rows)*nblk;
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
         zlacpy_("A", &rows_in_block, &cols_in_block, double_ptr, &na_rows, Buf_pos, &rows_in_block);     
         Buf_pos = Buf_pos + rows_in_block*cols_in_block;                         
         Size_send_U = Size_send_U + rows_in_block*cols_in_block; 
         cols_in_buffer_U_my_initial = cols_in_buffer_U_my_initial + cols_in_block; 
      }
      curr_col_loc = curr_col_loc + nblk;      
      rows_in_block = rows_in_block + ratio*nblk;
   }
   rows_in_buffer_U_my_initial = rows_in_block - ratio*nblk;    
   *Buf_pos = (double _Complex)cols_in_buffer_U_my_initial; 
   Buf_pos = Buf_pos + 1; 
   *Buf_pos = (double _Complex)rows_in_buffer_U_my_initial; 
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
            MPI_Sendrecv(Buf_to_send_U, (int) Size_send_U, ((MPI_Datatype)1275072546), (int) where_to_send_U, 0, Buf_U, (int) Size_U_stored, ((MPI_Datatype)1275072546), (int) from_where_to_receive_U, 0, row_comm, &status);
            MPI_Get_count(&status, ((MPI_Datatype)1275072546), &Size_receive_U_nowMPI);
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
            
         
         intNumber = ceil((double _Complex)cols_in_buffer_U_now/(double _Complex)nblk);  
         if(from_where_to_receive_U >= my_prow)
            rows_in_block = ceil(((double _Complex)(from_where_to_receive_U + 1) - (double _Complex)my_prow)/(double _Complex)np_rows)*nblk;  
         else
            rows_in_block = ratio*nblk; 
         for(j = 0; j < intNumber; j++)
         {
            if((j+1)*nblk < cols_in_buffer_U_now)
               cols_in_block = nblk; 
            else
               cols_in_block = cols_in_buffer_U_now - j*nblk;
               
            zlacpy_("A", &rows_in_block, &cols_in_block, CopyFrom, &rows_in_block, CopyTo, &rows_in_block);

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
            MPI_Sendrecv(Buf_to_send_U, (int) Size_send_U, ((MPI_Datatype)1275072546), (int) where_to_send_U, 0, Buf_to_receive_U, (int) Size_U_stored, ((MPI_Datatype)1275072546), (int) from_where_to_receive_U, 0, row_comm, &status);
            MPI_Get_count(&status, ((MPI_Datatype)1275072546), &Size_receive_UMPI);
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
         
         zlacpy_("A", &na_rows, &nb_cols, B, &na_rows, Buf_to_send_B, &na_rows);
         MPI_Sendrecv(Buf_to_send_B, (int) (nb_cols*na_rows), ((MPI_Datatype)1275072546), (int) where_to_send_B, 0, Buf_to_receive_B, (int) (nb_cols*Buf_rows), ((MPI_Datatype)1275072546), (int) from_where_to_receive_B, 0, col_comm, &status); 
         MPI_Get_count(&status, ((MPI_Datatype)1275072546), &Size_receive_BMPI); 
         Size_receive_B = (int) Size_receive_BMPI;
         Size_receive_B = Size_receive_B/nb_cols;    
	 
      }
      else
      {
         zlacpy_("A", &na_rows, &nb_cols, B, &na_rows, Buf_to_receive_B, &na_rows); 
         Size_receive_B = na_rows;
      }
   }
   else
   {
      zlacpy_("A", &na_rows, &nb_cols, B, &na_rows, Buf_to_receive_B, &na_rows);        
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
        
      
      MPI_Isend(Buf_to_send_U, (int) Size_send_U, ((MPI_Datatype)1275072546), (int) where_to_send_U, 0, row_comm, &request_U_Send); 
      MPI_Irecv(Buf_to_receive_U, (int) (ratio*Size_U_stored), ((MPI_Datatype)1275072546), (int) from_where_to_receive_U, 0, row_comm, &request_U_Recv);      
      
      MPI_Isend(Buf_to_send_B, (int) (Size_send_B*nb_cols), ((MPI_Datatype)1275072546), (int) where_to_send_B, 0, col_comm, &request_B_Send); 
      MPI_Irecv(Buf_to_receive_B, (int) (Buf_rows*nb_cols), ((MPI_Datatype)1275072546), (int) from_where_to_receive_B, 0, col_comm, &request_B_Recv);      
      
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
      
      num_of_blocks_in_U_buffer = ceil((double _Complex)cols_in_buffer_U/(double _Complex)nblk); 
      
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
      MPI_Get_count(&status, ((MPI_Datatype)1275072546), &Size_receive_UMPI); 
      Size_receive_U = (int) Size_receive_UMPI;

      MPI_Wait(&request_B_Send, &status);
      MPI_Wait(&request_B_Recv, &status);
      MPI_Get_count(&status, ((MPI_Datatype)1275072546), &Size_receive_BMPI); 
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
      
   num_of_blocks_in_U_buffer = ceil((double _Complex)cols_in_buffer_U/(double _Complex)nblk);
  
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


void cannons_triang_rectangular_c_dc(double _Complex* U, double _Complex* B, int local_rowsCast, int local_colsCast,
                                    int* u_desc, int* b_desc, double _Complex *Res, int row_comm, int col_comm)
{
  int local_rows, local_cols;

  local_rows = (int) local_rowsCast;
  local_cols = (int) local_colsCast;

  MPI_Comm c_row_comm = (MPI_Comm)(row_comm);
  MPI_Comm c_col_comm = (MPI_Comm)(col_comm);

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
















 
void cannons_reduction_c_dc(double _Complex* A, double _Complex* U, int local_rowsCast, int local_colsCasr, int* a_desc,
                            double _Complex *Res, int ToStore, int row_comm, int col_comm);















 
void cannons_triang_rectangular_c_dc(double _Complex* U, double _Complex* B, int local_rowsCast, int local_colsCast,
                                    int* u_desc, int* b_desc, double _Complex *Res, int row_comm, int col_comm);




















































































void dlacpy_(char*, int*, int*, double*, int*, double*, int*);
void dgemm_(char*, char*, int*, int*, int*, double*, double*, int*, double*, int*, double*, double*, int*); 


void slacpy_(char*, int*, int*, float*, int*, float*, int*);
void sgemm_(char*, char*, int*, int*, int*, float*, float*, int*, float*, int*, float*, float*, int*); 




void zlacpy_(char*, int*, int*, double _Complex*, int*, double _Complex*, int*);
void zgemm_(char*, char*, int*, int*, int*, double _Complex*, double _Complex*, int*, double _Complex*, int*, double _Complex*, double _Complex*, int*); 


void clacpy_(char*, int*, int*, float _Complex*, int*, float _Complex*, int*);
void cgemm_(char*, char*, int*, int*, int*, float _Complex*, float _Complex*, int*, float _Complex*, int*, float _Complex*, float _Complex*, int*); 




int numroc_(int*, int*, int*, int*, int*);


void pdlacpy_(char*, int*, int*, double*, int*, int*, int*, double*, int*, int*, int*);
void pdtran_(int*, int*, double*, double*, int*, int*, int*, double*, double*, int*, int*, int*);


void pslacpy_(char*, int*, int*, float*, int*, int*, int*, float*, int*, int*, int*);
void pstran_(int*, int*, float*, float*, int*, int*, int*, float*, float*, int*, int*, int*);



void pzlacpy_(char*, int*, int*, double _Complex*, int*, int*, int*, double _Complex*, int*, int*, int*);
void pztranc_(int*, int*, double _Complex*, double _Complex*, int*, int*, int*, double _Complex*, double _Complex*, int*, int*, int*);


void pclacpy_(char*, int*, int*, float _Complex*, int*, int*, int*, float _Complex*, int*, int*, int*);
void pctranc_(int*, int*, float _Complex*, float _Complex*, int*, int*, int*, float _Complex*, float _Complex*, int*, int*, int*);



void cannons_reduction_fc(float _Complex* A, float _Complex* U, int np_rows, int np_cols, int my_prow, int my_pcol,
                         int* a_desc, float _Complex *Res, int ToStore, MPI_Comm row_comm, MPI_Comm col_comm)
{
   
      
      
   
      
   
   
  
   int na, nblk, i, j, Size_send_A, Size_receive_A, Size_send_U, Size_receive_U, Buf_rows, Buf_cols, where_to_send_A, from_where_to_receive_A, where_to_send_U, from_where_to_receive_U, last_proc_row, last_proc_col, cols_in_buffer_A, rows_in_buffer_A, intNumber;
   float _Complex *Buf_to_send_A, *Buf_to_receive_A, *Buf_to_send_U, *Buf_to_receive_U, *data_ptr, *Buf_A, *Buf_pos, *U_local_start, *Res_ptr, *M, *M_T, *A_local_start, *U_local_start_curr, *U_stored, *CopyTo, *CopyFrom, *U_to_calc;
   int ratio, num_of_iters, cols_in_buffer, rows_in_block, rows_in_buffer, curr_col_loc, cols_in_block, curr_col_glob, curr_row_loc, Size_receive_A_now, Nb, owner, cols_in_buffer_A_now;
   int Size_receive_A_nowMPI, Size_receive_AMPI, Size_receive_UMPI;

   int  row_of_origin_U, rows_in_block_U, num_of_blocks_in_U_buffer, k, startPos, cols_in_buffer_U, rows_in_buffer_U, col_of_origin_A, curr_row_loc_res, curr_row_loc_A, curr_col_glob_res; 
   int curr_col_loc_res, curr_col_loc_buf, proc_row_curr, curr_col_loc_U, A_local_index, LDA_A, LDA_A_new, index_row_A_for_LDA, ii, rows_in_block_U_curr, width, row_origin_U, rows_in_block_A, cols_in_buffer_A_my_initial, rows_in_buffer_A_my_initial, proc_col_min;
   int *SizesU;
   int Size_U_skewed, Size_U_stored, Curr_pos_in_U_stored, rows_in_buffer_A_now;
   float _Complex done = 1.0;
   float _Complex dzero = 0.0;
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
   na_rows = numroc_(&na, &nblk, &my_prow, &zero, &np_rows);
   na_cols = numroc_(&na, &nblk, &my_pcol, &zero, &np_cols); 
   
   
   
   
   
   
   

   if (np_cols%np_rows != 0)
   {
      
      
      return;
   }
   if (np_cols < np_rows != 0)
   {
      
      
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
         Buf_rows = na_rows + 1;   
      else
         Buf_rows = na_rows + nblk;      
   else
      if(my_prow < last_proc_row)
         Buf_rows = na_rows;
      else if(my_prow > last_proc_row)
         Buf_rows = na_rows + nblk; 
      else  
         Buf_rows = na_rows + nblk - na_rows%nblk;  
      
   intNumber = ceil((float _Complex)na/(float _Complex)(np_cols*nblk));   
   Size_U_stored = ratio*nblk*nblk*intNumber*(intNumber+1)/2 + 2;   
   
   U_stored = malloc((Size_U_stored*(ToStore+1))*sizeof(float _Complex));
   SizesU = malloc(ToStore*sizeof(int));  
   Buf_to_send_A = malloc(ratio*Buf_cols*Buf_rows*sizeof(float _Complex));
   Buf_to_receive_A = malloc(ratio*Buf_cols*Buf_rows*sizeof(float _Complex));
   Buf_to_send_U = malloc(Size_U_stored*sizeof(float _Complex));
   Buf_to_receive_U = malloc(Size_U_stored*sizeof(float _Complex));
   if(ratio != 1)
      Buf_A = malloc(Buf_cols*Buf_rows*sizeof(float _Complex));   
   M = malloc(na_rows*na_cols*sizeof(float _Complex));
   M_T = malloc(na_rows*na_cols*sizeof(float _Complex));
   for(i = 0; i < na_rows*na_cols; i++)
      M[i] = 0; 
        
   
   
   
   if(ratio != 1)
      clacpy_("A", &na_rows, &na_cols, A, &na_rows, Buf_to_send_A, &na_rows);   
   Size_receive_A = 0; 
   
   
   for(i = 0; i < ratio; i++)
   {
      where_to_send_A = (my_pcol - my_prow - i*np_rows + np_cols)%np_cols;                
      from_where_to_receive_A = (my_pcol + my_prow + i*np_rows)%np_cols;
      
      
      if(ratio != 1)   
      {
         if(where_to_send_A != my_pcol)
         {
           MPI_Sendrecv(Buf_to_send_A, (int) (na_cols*na_rows), ((MPI_Datatype)1275070494),(int) where_to_send_A, (int) zero, Buf_A, (int) (na_rows*Buf_cols), ((MPI_Datatype)1275070494), (int) from_where_to_receive_A, (int) zero, row_comm, &status);
           MPI_Get_count(&status, ((MPI_Datatype)1275070494), &Size_receive_A_nowMPI);
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
         
         intNumber = ceil((float _Complex)Size_receive_A_now/(float _Complex)nblk);   
         for(j = 0; j < intNumber; j++)
         {
            width = nblk; 
            if(nblk*(j+1) > Size_receive_A_now)
               width = Size_receive_A_now - nblk*j; 
            clacpy_("A", &na_rows, &width, CopyFrom, &na_rows, CopyTo, &na_rows);
            CopyTo = CopyTo + na_rows*nblk*ratio; 
            CopyFrom = CopyFrom + na_rows*nblk; 
         }
      }
      else  
         if(my_prow > 0)
         {
            clacpy_("A", &na_rows, &na_cols, A, &na_rows, Buf_to_send_A, &na_rows);   
            MPI_Sendrecv(Buf_to_send_A, (int) (na_cols*na_rows), ((MPI_Datatype)1275070494), (int) where_to_send_A, (int) zero, Buf_to_receive_A, (int) (na_rows*Buf_cols), ((MPI_Datatype)1275070494), (int) from_where_to_receive_A, (int) zero, row_comm, &status);
            MPI_Get_count(&status, ((MPI_Datatype)1275070494), &Size_receive_AMPI);
            Size_receive_A = (int) Size_receive_AMPI;
	    Size_receive_A = Size_receive_A/na_rows;       
         }
         else
         {
            clacpy_("A", &na_rows, &na_cols, A, &na_rows, Buf_to_receive_A, &na_rows);   
            Size_receive_A = na_cols; 
         }
   }
   
   
     
   
   num_of_iters = ceil((float _Complex)na_cols/(float _Complex)nblk);             
   
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
      rows_in_block = ceil(((float _Complex)(my_pcol + 1) - (float _Complex)my_prow)/(float _Complex)np_rows)*nblk;
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
         clacpy_("A", &rows_in_block, &cols_in_block, data_ptr, &na_rows, Buf_pos, &rows_in_block);     
         Buf_pos = Buf_pos + rows_in_block*cols_in_block;                         
         Size_send_U = Size_send_U + rows_in_block*cols_in_block; 
      }
      curr_col_loc = curr_col_loc + nblk;      
      rows_in_block = rows_in_block + ratio*nblk;
   }
   rows_in_buffer = rows_in_block - ratio*nblk;    
   *Buf_pos = (float _Complex)rows_in_buffer; 
   Size_send_U = Size_send_U + 1;
   
   
   if(where_to_send_U != my_prow)
   {   
      
      MPI_Sendrecv(Buf_to_send_U, (int) Size_send_U, ((MPI_Datatype)1275070494), (int) where_to_send_U, (int) zero, Buf_to_receive_U, (int) (Buf_rows*na_cols), ((MPI_Datatype)1275070494), (int) from_where_to_receive_U, (int) zero, col_comm, &status); 
      MPI_Get_count(&status, ((MPI_Datatype)1275070494), &Size_receive_UMPI); 
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
      MPI_Isend(Buf_to_send_A, (int) (Size_send_A*na_rows), ((MPI_Datatype)1275070494), (int) where_to_send_A, (int) zero, row_comm, &request_A_Send); 
      MPI_Irecv(Buf_to_receive_A, (int) (Buf_cols*na_rows*ratio), ((MPI_Datatype)1275070494), (int) from_where_to_receive_A, (int) zero, row_comm, &request_A_Recv);
         
      
      Size_send_U = Size_receive_U; 
      MPI_Isend(Buf_to_send_U, (int) Size_send_U, ((MPI_Datatype)1275070494), (int) where_to_send_U, (int) zero, col_comm, &request_U_Send); 
      MPI_Irecv(Buf_to_receive_U, (int) (Buf_rows*na_cols), ((MPI_Datatype)1275070494), (int) from_where_to_receive_U, (int) zero, col_comm, &request_U_Recv); 
      
      
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
    
      num_of_blocks_in_U_buffer = ceil(((float _Complex)cols_in_buffer - (float _Complex)curr_col_loc_buf)/(float _Complex)nblk); 
      
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

      MPI_Get_count(&status, ((MPI_Datatype)1275070494), &Size_receive_AMPI); 
      Size_receive_A = (int) Size_receive_AMPI;
      Size_receive_A = Size_receive_A / na_rows;
      
      
      MPI_Wait(&request_U_Send, &status);
      MPI_Wait(&request_U_Recv, &status);
      MPI_Get_count(&status, ((MPI_Datatype)1275070494), &Size_receive_UMPI); 
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
    
   num_of_blocks_in_U_buffer = ceil(((float _Complex)cols_in_buffer - (float _Complex)curr_col_loc_buf)/(float _Complex)nblk); 
      
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
   
   
   
   pctranc_(&na, &na, &done, M, &one, &one, a_desc, &dzero, M_T, &one, &one, a_desc);     
 
   
           
   
   
   
   if((ratio != 1)||(my_prow != 0))   
      Buf_pos = Buf_to_send_A;     
   else
      Buf_pos = Buf_to_receive_A;  
   
   
   num_of_iters = ceil((float _Complex)na_cols/(float _Complex)nblk);             
   
   cols_in_buffer_A_my_initial = 0;
   Size_send_A = 0; 
   
   if(my_pcol <= my_prow)  
   {
      curr_row_loc = 0;     
      rows_in_buffer_A_my_initial = na_rows;
   }
   else
   {
      curr_row_loc = ceil((float _Complex)(((float _Complex)my_pcol - (float _Complex)my_prow)/(float _Complex)np_rows))*nblk; 
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
         clacpy_("A", &rows_in_block, &cols_in_block, A_local_start, &na_rows, Buf_pos, &rows_in_block);     
         Buf_pos = Buf_pos + rows_in_block*cols_in_block;
         Size_send_A = Size_send_A + rows_in_block*cols_in_block; 
         cols_in_buffer_A_my_initial = cols_in_buffer_A_my_initial + cols_in_block; 
      }
      curr_row_loc = curr_row_loc + ratio*nblk;
   }
   *Buf_pos = (float _Complex)cols_in_buffer_A_my_initial; 
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
            MPI_Sendrecv(Buf_to_send_A, (int) Size_send_A, ((MPI_Datatype)1275070494), (int) where_to_send_A, (int) zero, Buf_A, (int) Size_U_stored, ((MPI_Datatype)1275070494), (int) from_where_to_receive_A, (int) zero, row_comm, &status);
            MPI_Get_count(&status, ((MPI_Datatype)1275070494), &Size_receive_A_nowMPI);
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
               rows_in_buffer_A_now = na_rows - ceil((float _Complex)(((float _Complex)from_where_to_receive_A - (float _Complex)my_prow)/(float _Complex)np_rows))*nblk; 
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
            
         
         intNumber = ceil((float _Complex)cols_in_buffer_A_now/(float _Complex)nblk);  
         rows_in_block = rows_in_buffer_A_now; 
         for(j = 0; j < intNumber; j++)
         {
            if((j+1)*nblk < cols_in_buffer_A_now)
               cols_in_block = nblk; 
            else
               cols_in_block = cols_in_buffer_A_now - j*nblk;
               
            clacpy_("A", &rows_in_block, &cols_in_block, CopyFrom, &rows_in_block, CopyTo, &rows_in_block);

            CopyFrom = CopyFrom + rows_in_block*cols_in_block; 
            CopyTo = CopyTo + nblk*(ratio*rows_in_block - nblk*(ratio-1)*ratio/2);  
            rows_in_block = rows_in_block - ratio*nblk;     
         }
      }
      else    
      {
         if(my_prow > 0)
         {
            MPI_Sendrecv(Buf_to_send_A, (int) Size_send_A, ((MPI_Datatype)1275070494), (int) where_to_send_A, (int) zero, Buf_to_receive_A, (int) Size_U_stored, ((MPI_Datatype)1275070494), (int) from_where_to_receive_A, (int) zero, row_comm, &status);
            MPI_Get_count(&status, ((MPI_Datatype)1275070494), &Size_receive_AMPI);
            Size_receive_A = (int) Size_receive_AMPI;

            cols_in_buffer_A = (int)Buf_to_receive_A[Size_receive_A-1];
            if(from_where_to_receive_A <= my_prow)  
            {
               rows_in_buffer_A = na_rows;
            }
            else
            {
               rows_in_buffer_A = na_rows - ceil((float _Complex)(((float _Complex)from_where_to_receive_A - (float _Complex)my_prow)/(float _Complex)np_rows))*nblk; 
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
      MPI_Isend(Buf_to_send_A, (int) Size_send_A, ((MPI_Datatype)1275070494), (int) where_to_send_A, (int) zero, row_comm, &request_A_Send); 
      MPI_Irecv(Buf_to_receive_A, (int) (ratio*Size_U_stored), ((MPI_Datatype)1275070494), (int) from_where_to_receive_A, (int) zero, row_comm, &request_A_Recv);
         
      
      Size_send_U = Size_receive_U; 
      if (j > ToStore)
      {
         if(j > ToStore + 1)
         {
            MPI_Isend(Buf_to_send_U, (int) Size_send_U, ((MPI_Datatype)1275070494), (int) where_to_send_U, (int) zero, col_comm, &request_U_Send);
            U_to_calc = Buf_to_send_U;
         }
         else {
	    MPI_Isend(U_to_calc, (int) Size_send_U, ((MPI_Datatype)1275070494), (int) where_to_send_U, (int) zero, col_comm, &request_U_Send);
	 }
         MPI_Irecv(Buf_to_receive_U, (int) Size_U_stored, ((MPI_Datatype)1275070494), (int) from_where_to_receive_U, (int) zero, col_comm, &request_U_Recv);	 
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
      
      num_of_blocks_in_U_buffer = ceil((float _Complex)((float _Complex)cols_in_buffer_U/(float _Complex)nblk)); 
      if(my_pcol >= row_of_origin_U)    
         rows_in_block_U = ceil(((float _Complex)(my_pcol + 1) - (float _Complex)row_of_origin_U)/(float _Complex)np_rows)*nblk;  
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
 
            
            for (ii = 0; ii < ceil((float _Complex)rows_in_block_U/(float _Complex)nblk); ii++)
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
      MPI_Get_count(&status, ((MPI_Datatype)1275070494), &Size_receive_AMPI); 
      Size_receive_A = (int) Size_receive_AMPI;
      
      if (j <= ToStore)
      {
         U_to_calc = &U_stored[Curr_pos_in_U_stored];
         Curr_pos_in_U_stored = Curr_pos_in_U_stored + SizesU[j-1]; 
         Size_receive_U =  SizesU[j-1];
      }
      else
      {
         MPI_Wait(&request_U_Send, &status);
         MPI_Wait(&request_U_Recv, &status);
	 MPI_Get_count(&status, ((MPI_Datatype)1275070494), &Size_receive_UMPI); 
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
      
   num_of_blocks_in_U_buffer = ceil((float _Complex)((float _Complex)cols_in_buffer_U/(float _Complex)nblk));
   if(my_pcol >= row_of_origin_U)    
      rows_in_block_U = ceil(((float _Complex)(my_pcol + 1) - (float _Complex)row_of_origin_U)/(float _Complex)np_rows)*nblk;  
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

         
         for (ii = 0; ii < ceil((float _Complex)rows_in_block_U/(float _Complex)nblk); ii++)
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
   
   pctranc_(&na, &na, &done, Res, &one, &one, a_desc, &dzero, M, &one, &one, a_desc);
   pclacpy_("U", &na, &na, M, &one, &one, a_desc, Res, &one, &one, a_desc);
      

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

void cannons_reduction_c_fc(float _Complex* A, float _Complex* U, int local_rowsCast, int local_colsCast,
                         int* a_desc, float _Complex *Res, int ToStore, int row_comm, int col_comm)
{
  int local_rows, local_cols;
  local_rows = (int) local_rowsCast;
  local_cols = (int) local_colsCast;

  MPI_Comm c_row_comm = (MPI_Comm)(row_comm);
  MPI_Comm c_col_comm = (MPI_Comm)(col_comm);


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





























































void cannons_triang_rectangular_fc(float _Complex* U, float _Complex* B, int np_rows, int np_cols, int my_prow, int my_pcol, int* U_desc, int* b_desc, float _Complex *Res, MPI_Comm row_comm, MPI_Comm col_comm)
{
   
   
   
   
   
   
   
   
   
  
   int na, nb, nblk, width, na_rows, na_cols, nb_cols, cols_in_buffer_U_my_initial, cols_in_buffer_U, rows_in_buffer_U, Size_receive_U_now, rows_in_buffer_U_now, cols_in_buffer_U_now, rows_in_buffer_U_my_initial;

   int Size_receive_U_nowMPI, Size_receive_UMPI, Size_receive_BMPI;
   int i, j, Size_send_U, Size_receive_U, Size_send_B, Size_receive_B, intNumber, Buf_rows, Buf_cols_U, Buf_cols_B, curr_rows, num_of_iters, cols_in_buffer, rows_in_block, curr_col_loc, cols_in_block, num_of_blocks_in_U_buffer, col_of_origin_U, b_rows_mult, b_cols_mult; 
   
   float _Complex *Buf_to_send_U, *Buf_to_receive_U, *Buf_to_send_B, *Buf_to_receive_B, *Buf_U, *PosBuff;
  
   int where_to_send_U, from_where_to_receive_U, where_to_send_B, from_where_to_receive_B, last_proc_col_B, last_proc_row_B, n, Size_U_stored, proc_col_min; 
   
   float _Complex *U_local_start, *Buf_pos, *B_local_start, *double_ptr, *CopyTo, *CopyFrom;
   
   int ratio;
   
   MPI_Status status;

   int one = 1;
   int zero = 0; 
   float _Complex done = 1.0;
   float _Complex dzero = 0.0;
      
   na = U_desc[2];
   nblk = U_desc[4]; 
   nb = b_desc[3];
   
   na_rows = numroc_(&na, &nblk, &my_prow, &zero, &np_rows);
   na_cols = numroc_(&na, &nblk, &my_pcol, &zero, &np_cols);
   nb_cols = numroc_(&nb, &nblk, &my_pcol, &zero, &np_cols);
   
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
   
   intNumber = ceil((float _Complex)na/(float _Complex)(np_cols*nblk));   
   Size_U_stored = ratio*nblk*nblk*intNumber*(intNumber+1)/2 + 2;   
   
   Buf_to_send_U = malloc(ratio*Size_U_stored*sizeof(float _Complex));
   Buf_to_receive_U = malloc(ratio*Size_U_stored*sizeof(float _Complex));
   Buf_to_send_B = malloc(Buf_cols_B*Buf_rows*sizeof(float _Complex));
   Buf_to_receive_B = malloc(Buf_cols_B*Buf_rows*sizeof(float _Complex));
   if(ratio != 1)
      Buf_U = malloc(Size_U_stored*sizeof(float _Complex));   
    
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
      
   num_of_iters = ceil((float _Complex)na_cols/(float _Complex)nblk);             
   num_of_iters = num_of_iters - curr_col_loc;   
   curr_col_loc = curr_col_loc*nblk;             

   if(my_pcol >= my_prow )
      rows_in_block = ceil(((float _Complex)(my_pcol + 1) - (float _Complex)my_prow)/(float _Complex)np_rows)*nblk;
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
         clacpy_("A", &rows_in_block, &cols_in_block, double_ptr, &na_rows, Buf_pos, &rows_in_block);     
         Buf_pos = Buf_pos + rows_in_block*cols_in_block;                         
         Size_send_U = Size_send_U + rows_in_block*cols_in_block; 
         cols_in_buffer_U_my_initial = cols_in_buffer_U_my_initial + cols_in_block; 
      }
      curr_col_loc = curr_col_loc + nblk;      
      rows_in_block = rows_in_block + ratio*nblk;
   }
   rows_in_buffer_U_my_initial = rows_in_block - ratio*nblk;    
   *Buf_pos = (float _Complex)cols_in_buffer_U_my_initial; 
   Buf_pos = Buf_pos + 1; 
   *Buf_pos = (float _Complex)rows_in_buffer_U_my_initial; 
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
            MPI_Sendrecv(Buf_to_send_U, (int) Size_send_U, ((MPI_Datatype)1275070494), (int) where_to_send_U, 0, Buf_U, (int) Size_U_stored, ((MPI_Datatype)1275070494), (int) from_where_to_receive_U, 0, row_comm, &status);
            MPI_Get_count(&status, ((MPI_Datatype)1275070494), &Size_receive_U_nowMPI);
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
            
         
         intNumber = ceil((float _Complex)cols_in_buffer_U_now/(float _Complex)nblk);  
         if(from_where_to_receive_U >= my_prow)
            rows_in_block = ceil(((float _Complex)(from_where_to_receive_U + 1) - (float _Complex)my_prow)/(float _Complex)np_rows)*nblk;  
         else
            rows_in_block = ratio*nblk; 
         for(j = 0; j < intNumber; j++)
         {
            if((j+1)*nblk < cols_in_buffer_U_now)
               cols_in_block = nblk; 
            else
               cols_in_block = cols_in_buffer_U_now - j*nblk;
               
            clacpy_("A", &rows_in_block, &cols_in_block, CopyFrom, &rows_in_block, CopyTo, &rows_in_block);

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
            MPI_Sendrecv(Buf_to_send_U, (int) Size_send_U, ((MPI_Datatype)1275070494), (int) where_to_send_U, 0, Buf_to_receive_U, (int) Size_U_stored, ((MPI_Datatype)1275070494), (int) from_where_to_receive_U, 0, row_comm, &status);
            MPI_Get_count(&status, ((MPI_Datatype)1275070494), &Size_receive_UMPI);
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
         
         clacpy_("A", &na_rows, &nb_cols, B, &na_rows, Buf_to_send_B, &na_rows);
         MPI_Sendrecv(Buf_to_send_B, (int) (nb_cols*na_rows), ((MPI_Datatype)1275070494), (int) where_to_send_B, 0, Buf_to_receive_B, (int) (nb_cols*Buf_rows), ((MPI_Datatype)1275070494), (int) from_where_to_receive_B, 0, col_comm, &status); 
         MPI_Get_count(&status, ((MPI_Datatype)1275070494), &Size_receive_BMPI); 
         Size_receive_B = (int) Size_receive_BMPI;
         Size_receive_B = Size_receive_B/nb_cols;    
	 
      }
      else
      {
         clacpy_("A", &na_rows, &nb_cols, B, &na_rows, Buf_to_receive_B, &na_rows); 
         Size_receive_B = na_rows;
      }
   }
   else
   {
      clacpy_("A", &na_rows, &nb_cols, B, &na_rows, Buf_to_receive_B, &na_rows);        
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
        
      
      MPI_Isend(Buf_to_send_U, (int) Size_send_U, ((MPI_Datatype)1275070494), (int) where_to_send_U, 0, row_comm, &request_U_Send); 
      MPI_Irecv(Buf_to_receive_U, (int) (ratio*Size_U_stored), ((MPI_Datatype)1275070494), (int) from_where_to_receive_U, 0, row_comm, &request_U_Recv);      
      
      MPI_Isend(Buf_to_send_B, (int) (Size_send_B*nb_cols), ((MPI_Datatype)1275070494), (int) where_to_send_B, 0, col_comm, &request_B_Send); 
      MPI_Irecv(Buf_to_receive_B, (int) (Buf_rows*nb_cols), ((MPI_Datatype)1275070494), (int) from_where_to_receive_B, 0, col_comm, &request_B_Recv);      
      
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
      
      num_of_blocks_in_U_buffer = ceil((float _Complex)cols_in_buffer_U/(float _Complex)nblk); 
      
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
      MPI_Get_count(&status, ((MPI_Datatype)1275070494), &Size_receive_UMPI); 
      Size_receive_U = (int) Size_receive_UMPI;

      MPI_Wait(&request_B_Send, &status);
      MPI_Wait(&request_B_Recv, &status);
      MPI_Get_count(&status, ((MPI_Datatype)1275070494), &Size_receive_BMPI); 
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
      
   num_of_blocks_in_U_buffer = ceil((float _Complex)cols_in_buffer_U/(float _Complex)nblk);
  
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


void cannons_triang_rectangular_c_fc(float _Complex* U, float _Complex* B, int local_rowsCast, int local_colsCast,
                                    int* u_desc, int* b_desc, float _Complex *Res, int row_comm, int col_comm)
{
  int local_rows, local_cols;

  local_rows = (int) local_rowsCast;
  local_cols = (int) local_colsCast;

  MPI_Comm c_row_comm = (MPI_Comm)(row_comm);
  MPI_Comm c_col_comm = (MPI_Comm)(col_comm);

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

















 

void cannons_reduction_c_fc(float _Complex* A, float _Complex* U, int local_rowsCast, int local_colsCast, int* a_desc,
                         float _Complex *Res, int ToStore, int row_comm, int col_comm);















 
void cannons_triang_rectangular_c_fc(float _Complex* U, float _Complex* B, int local_rowsCast, int local_colsCast,
                                    int* u_desc, int* b_desc, float _Complex *Res, int row_comm, int col_comm);
