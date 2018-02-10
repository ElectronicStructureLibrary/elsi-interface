#!/usr/bin/perl

# Copyright (c) 2015-2018, the ELSI team. All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#
#  * Redistributions of source code must retain the above copyright notice,
#    this list of conditions and the following disclaimer.
#
#  * Redistributions in binary form must reproduce the above copyright notice,
#    this list of conditions and the following disclaimer in the documentation
#    and/or other materials provided with the distribution.
#
#  * Neither the name of the "ELectronic Structure Infrastructure" project nor
#    the names of its contributors may be used to endorse or promote products
#    derived from this software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
# ARE DISCLAIMED. IN NO EVENT SHALL COPYRIGHT HOLDER BE LIABLE FOR ANY DIRECT,
# INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
# BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
# DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY
# OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
# NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE,
# EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

# This script is based on a similar script from FHI-aims, which was modified to
# better fit into the ELSI design philosophy.

use strict;

# script to write a Fortran subroutine which writes the Version stamp and
# compile time.

# Obtain version stamp from input arg:
my $RELEASE_DATE = $ARGV[0];

# obtain local time and date from PERL function, and fix up to get proper format
my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime(time);
$year = $year + 1900;
$mon  = $mon + 1;
my $date = sprintf("%4d-%02d-%02d", $year, $mon, $mday);
my $time = sprintf("%02d:%02d:%02d", $hour, $min, $sec);

# now do the same thing, but in UTC and generate a RFC3339-formatted timestamp
my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = gmtime(time);
$year = $year + 1900;
$mon  = $mon + 1;
my $datetime = sprintf("%4d-%02d-%02dT%02d:%02d:%02dZ", $year, $mon, $mday, $hour, $min, $sec);

# Where am I?
my $hostname = qx/hostname/;
chomp $hostname;
my $host = "$hostname";

# Can we get the current Git ID?
my ($isgitrepos, $gitmsg, $gitret, $ismod, $modlog, $modstr, $gitcommit, $gitcommitabbrev, $gitrev, $gitstr);

# Determine whether the git commit has been modified
my $perlret = system("git diff --raw --exit-code HEAD >/dev/null 2>&1");
if ($perlret == -1) {
    # execution failed quite early
    $isgitrepos = 0;
} else {
    print STDERR $gitret;
    $gitret = $perlret >> 8;
    $isgitrepos = $gitret == 0 || $gitret == 1;
}

if ($isgitrepos) {
    # Put the status of the git commit's integrity into pretty strings
    $ismod = $gitret == 1;
    if ($ismod) { $modstr = " (modified)"; } else { $modstr = ""; }
    if ($ismod) { $modlog = ".true."; } else { $modlog = ".false."; }

    # WPH: Get the git commit
    $gitrev = qx/git show --no-color --pretty=oneline | head -n 1/;
    chomp $gitrev;
    # CC: clean git string from special chars that can break the compilation
    $gitrev =~ s/\@/[at]/g;
    $gitrev =~ s/[^A-Za-z 0-9\.,\[\]:\/-_]//g;
    if ( length($gitrev) > 40){
    	$gitcommit = substr( $gitrev, 0, 40);
    	$gitcommitabbrev = substr( $gitrev, 0, 7);
    };

    # WPH: Get only the git commit message
    $gitrev = qx/git show --no-color --pretty=oneline --abbrev-commit | head -n 1/;
    chomp $gitrev;
    # CC: clean git string from special chars that can break the compilation
    $gitrev =~ s/\@/[at]/g;
    $gitrev =~ s/[^A-Za-z 0-9\.,\[\]:\/-_]//g;
    if ( length($gitrev) > 48){
    	$gitrev = substr( $gitrev, 0, 48);
    };
    $gitmsg = substr( $gitrev, 8, length($gitrev));

# WPH: The following is the original functionality for outputting a git commit + message.
#      I've kept it here in case we ever want to use a similar format, and because
#      distributing a git-rev.txt file could be useful to get around the requirement that
#      users have git installed, so I figured I'd leave the code in just in case.
#    $gitrev = qx/git show --no-color --pretty=oneline --abbrev-commit | head -n 1/;
#    chomp $gitrev;
#    # CC: clean git string from special chars that can break the compilation
#    $gitrev =~ s/\@/[at]/g;
#    $gitrev =~ s/[^A-Za-z 0-9\.,\[\]:\/-_]//g;
#    if ( length($gitrev) > 50){
#    	$gitrev = substr( $gitrev, 0, 50);
#	$gitrev = $gitrev."[...]";
#    };
#    $gitstr = "Git rev.$modstr: $gitrev";
#} elsif (-f "git-rev.txt") {
#    my $catgitstr = qx/cat git-rev.txt/;
#    chomp $catgitstr;
#    $gitstr = "Based on $catgitstr";
} else {
    $gitcommit = "UNKNOWN";
    $gitcommitabbrev = "UNKNOWN";
    $modlog = ".false.";
    $gitmsg = "UNKNOWN";
}

if ($RELEASE_DATE eq "--only-gitstr") {
    print "$gitstr\n";
} else {
    print "! Copyright (c) 2015-2018, the ELSI team. All rights reserved.\n";
    print "!\n";
    print "! Redistribution and use in source and binary forms, with or without\n";
    print "! modification, are permitted provided that the following conditions are met:\n";
    print "!\n";
    print "!  * Redistributions of source code must retain the above copyright notice,\n";
    print "!    this list of conditions and the following disclaimer.\n";
    print "!\n";
    print "!  * Redistributions in binary form must reproduce the above copyright notice,\n";
    print "!    this list of conditions and the following disclaimer in the documentation\n";
    print "!    and/or other materials provided with the distribution.\n";
    print "!\n";
    print "!  * Neither the name of the \"ELectronic Structure Infrastructure\" project nor\n";
    print "!    the names of its contributors may be used to endorse or promote products\n";
    print "!    derived from this software without specific prior written permission.\n";
    print "!\n";
    print "! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS \"AS IS\"\n";
    print "! AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE\n";
    print "! IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE\n";
    print "! ARE DISCLAIMED. IN NO EVENT SHALL COPYRIGHT HOLDER BE LIABLE FOR ANY DIRECT,\n";
    print "! INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,\n";
    print "! BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,\n";
    print "! DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY\n";
    print "! OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING\n";
    print "! NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE,\n";
    print "! EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.\n";
    print "\n";
    print "!>\n";
    print "!! This module provides details about ELSI's compilation, generated by the\n";
    print "!! script elsi_version_module_writer.pl as part of the ELSI make process.\n";
    print "!!\n";
    print "module ELSI_VERSION\n";
    print "\n";
    print "  implicit none\n";
    print "\n";
    print "  public\n";
    print "\n";
    print "  character(len=10), parameter :: RELEASE_DATE            = \"$RELEASE_DATE\"\n";
    print "  character(len=40), parameter :: GIT_COMMIT              = \"$gitcommit\"\n";
    print "  character(len=7),  parameter :: GIT_COMMIT_ABBREV       = \"$gitcommitabbrev\"\n";
    print "  logical,           parameter :: GIT_COMMIT_WAS_MODIFIED = $modlog\n";
    print "  character(len=40), parameter :: GIT_COMMIT_MSG_ABBREV   = \"$gitmsg\"\n";
    print "  character(len=*),  parameter :: SOURCE_HOSTNAME         = \"$host\"\n";
    print "  character(len=10), parameter :: SOURCE_LOCAL_DATE       = \"$date\"\n";
    print "  character(len=8),  parameter :: SOURCE_LOCAL_TIME       = \"$time\"\n";
    print "  character(len=20), parameter :: SOURCE_DATETIME         = \"$datetime\"\n";
    print "\n";
    print "end module ELSI_VERSION\n";
}
