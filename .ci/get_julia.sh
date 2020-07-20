#!/bin/bash
set -eu

LINUX_URL=https://julialang-s3.julialang.org/bin/linux/x64/1.4/julia-1.4.2-linux-x86_64.tar.gz
LINUX_HASH=d77311be23260710e89700d0b1113eecf421d6cf31a9cebad3f6bdd606165c28

if [ "$TRAVIS_OS_NAME" == "osx" ]; then
	:  # Done in .travis.yml
	# export HOMEBREW_NO_INSTALL_CLEANUP=1
	# export HOMEBREW_NO_AUTO_UPDATE=1
	# brew cask install julia
else
	mkdir -p "$HOME/julia_binary"
	pushd "$HOME/julia_binary"

	wget "$LINUX_URL" -O julia.tar.gz
	SHASUM=$(sha256sum julia.tar.gz | cut -f1 -d" ")
	if [ "$SHASUM" != "$LINUX_HASH" ]; then
		echo "Julia hashsum mismatch" >&2
		exit 1
	fi

	tar xzf julia.tar.gz
	JULIAPATH=$(echo $PWD/julia-*/bin)/julia
	echo "Installed julia to $JULIAPATH"
	popd

	pushd "$HOME/bin"
	ln -s $JULIAPATH julia
	popd
fi

