#!/bin/bash
set -eu

LINUX_URL=https://julialang-s3.julialang.org/bin/linux/x64/1.4/julia-1.4.0-linux-x86_64.tar.gz
LINUX_HASH=30d126dc3598f3cd0942de21cc38493658037ccc40eb0882b3b4c418770ca751

if [ "$TRAVIS_OS_NAME" == "osx" ]; then
	export HOMEBREW_NO_INSTALL_CLEANUP=1
	export HOMEBREW_NO_AUTO_UPDATE=1
	brew cask install julia
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

